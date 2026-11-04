"""
FMUMixin with the FMU of examples/basic (OpenModelica, fmi-ls-dae). The FMU is a build artifact
that is not in the repository: export it with OpenModelica (``--daeMode --simCodeTarget=C+Rust
--fmiInputDerivatives``, ``buildModelFMU(Example, version="3.0", fmuType="me")``) to
``examples/basic/model/Example.fmu``, or point ``RTCTOOLS_TEST_FMU`` at it. The tests are skipped
without it, and without a CasADi with fmi-ls-dae support.
"""

import os
import unittest

import casadi as ca
import numpy as np

from rtctools.optimization.collocated_integrated_optimization_problem import (
    CollocatedIntegratedOptimizationProblem,
)
from rtctools.optimization.csv_mixin import CSVMixin

from ..test_case import TestCase

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "..", "examples", "basic")
FMU = os.environ.get("RTCTOOLS_TEST_FMU", os.path.join(EXAMPLE, "model", "Example.fmu"))


def _skip_reason():
    if not os.path.exists(FMU):
        return f"FMU {FMU} not found"
    if not hasattr(ca, "DaeBuilder"):
        return "CasADi without DaeBuilder"
    return None


@unittest.skipIf(_skip_reason() is not None, _skip_reason() or "")
class TestFMUMixin(TestCase):
    @classmethod
    def setUpClass(cls):
        from rtctools.optimization.fmu_mixin import FMUMixin

        class Example(CSVMixin, FMUMixin, CollocatedIntegratedOptimizationProblem):
            model_name = "Example"

            def objective(self, ensemble_member):
                return self.integral("Q_release", ensemble_member=ensemble_member)

        cls.problem = Example(
            model_folder=os.path.dirname(FMU),
            input_folder=os.path.join(EXAMPLE, "input"),
            output_folder=os.path.join(EXAMPLE, "output"),
        )

    def names(self, key):
        return [v.name() for v in self.problem.dae_variables[key]]

    def test_variables(self):
        self.assertEqual(self.names("states"), ["storage.V"])
        self.assertEqual(self.names("derivatives"), ["der(storage.V)"])
        self.assertIn("V_storage", self.names("algebraics"))
        self.assertEqual(self.names("control_inputs"), ["Q_release"])
        self.assertEqual(self.names("constant_inputs"), ["Q_in"])
        self.assertEqual(self.names("time"), ["time"])
        # FMI <Alias> elements
        self.assertEqual(
            self.problem.alias_relation.canonical_signed("storage.Q_release"), ("Q_release", 1)
        )
        # Bounds and nominals from the FMU attributes
        bounds = self.problem.bounds()
        self.assertEqual(bounds["storage.V"], (2e5, 6e5))
        self.assertEqual(bounds["Q_release"], (0.0, 6.5))
        self.assertEqual(self.problem.variable_nominal("storage.V"), 4e5)
        # Parameters (start values)
        parameters = self.problem.parameters(0)
        self.assertEqual(parameters["storage.Q_nominal"], 1.0)

    def test_jacobian_directional_derivatives(self):
        """The residual Jacobian (FMI directional derivatives) matches finite differences."""
        dv = self.problem.dae_variables
        syms = dv["states"] + dv["derivatives"] + dv["algebraics"] + dv["control_inputs"]
        syms += dv["constant_inputs"]
        x = ca.vertcat(*syms)
        J = ca.Function("J", [dv["time"][0], x], [ca.jacobian(self.problem.dae_residual, x)])
        f = ca.Function("f", [dv["time"][0], x], [self.problem.dae_residual])
        rng = np.random.default_rng(3)
        big = ("storage.V", "V_storage")
        x0 = np.array([4e5 if s.name() in big else rng.uniform(1, 3) for s in syms])
        Ja = np.array(J(0.0, x0))
        Jd = np.zeros_like(Ja)
        for j in range(x0.size):
            h = 1e-6 * max(1.0, abs(x0[j]))
            e = np.zeros_like(x0)
            e[j] = h
            Jd[:, j] = (np.array(f(0.0, x0 + e)) - np.array(f(0.0, x0 - e))).ravel() / (2 * h)
        np.testing.assert_allclose(Ja, Jd, rtol=1e-6, atol=1e-6)

    def test_optimize(self):
        self.problem.optimize()
        self.assertTrue(self.problem.solver_stats["success"])
        results = self.problem.extract_results()
        np.testing.assert_allclose(results["V_storage"], results["storage.V"], rtol=1e-8)

    def test_residual_against_pymoca(self):
        """
        Parity with ModelicaMixin (pymoca) on the same Modelica model: at random points, matched
        by variable name (and aliases), each FMU residual row is a multiple of a pymoca row.
        """
        from rtctools.optimization.modelica_mixin import ModelicaMixin

        class ExamplePymoca(CSVMixin, ModelicaMixin, CollocatedIntegratedOptimizationProblem):
            model_name = "Example"

        try:
            pm = ExamplePymoca(
                model_folder=os.path.join(EXAMPLE, "model"),
                input_folder=os.path.join(EXAMPLE, "input"),
                output_folder=os.path.join(EXAMPLE, "output"),
            )
        except Exception as e:  # e.g. Modelica Standard Library not on MODELICAPATH
            self.skipTest(f"pymoca model not available: {e}")

        def residual_function(problem):
            dv = problem.dae_variables
            syms = dv["time"] + dv["states"] + dv["derivatives"] + dv["algebraics"]
            syms += dv["control_inputs"] + dv["constant_inputs"]
            p = problem.parameters(0)
            res = ca.substitute(
                problem.dae_residual,
                ca.vertcat(*dv["parameters"]),
                ca.vertcat(*[float(p[s.name()]) for s in dv["parameters"]]),
            )
            return ca.Function("f", syms, [res]), [s.name() for s in syms]

        f_pm, n_pm = residual_function(pm)
        f_fmu, n_fmu = residual_function(self.problem)
        rng = np.random.default_rng(2)
        rows_pm, rows_fmu = [], []
        for _ in range(10):
            point = {n: rng.uniform(-10, 10) for n in n_pm}
            point["storage.V"] = rng.uniform(2e5, 6e5)

            def value(n, point=point):
                canonical, sign = pm.alias_relation.canonical_signed(n)
                return sign * point[canonical]

            rows_pm.append(np.array(f_pm(*[point[n] for n in n_pm])).ravel())
            rows_fmu.append(np.array(f_fmu(*[value(n) for n in n_fmu])).ravel())
        rows_pm, rows_fmu = np.array(rows_pm), np.array(rows_fmu)
        for i in range(rows_fmu.shape[1]):
            if np.allclose(rows_fmu[:, i], 0.0, atol=1e-8):
                continue  # e.g. V_storage - F_y(storage.V), with V_storage aliased by pymoca
            c, *_ = np.linalg.lstsq(rows_pm, rows_fmu[:, i], rcond=None)
            np.testing.assert_allclose(rows_pm @ c, rows_fmu[:, i], rtol=1e-12, atol=1e-9)
