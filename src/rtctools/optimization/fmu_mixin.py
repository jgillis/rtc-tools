import atexit
import csv
import logging
import os
import re

import casadi as ca
import numpy as np

from rtctools._internal.alias_tools import AliasDict, AliasRelation
from rtctools._internal.caching import cached
from rtctools._internal.ensemble_bounds_decorator import ensemble_bounds_check
from rtctools._internal.fmu_description import FMUDescription

from .optimization_problem import OptimizationProblem
from .timeseries import Timeseries

logger = logging.getLogger("rtctools")

# FMU-internal variables that, like in ModelicaMixin's default pymoca options, are not exposed
_ELIMINABLE = re.compile(r"(.*[.]|^)_\w+(\[[\d,]+\])?\Z")

# Options of the FmuFunctions for the DAE residuals: forward directional derivatives of the FMU
# (fmi3GetDirectionalDerivative); its adjoint derivatives are not used (see FMUMixin docstring).
_NO_ADJOINT = {"uses_adjoint_derivatives": False}
# Options of the FmuFunction for FMU outputs: finite differences, as the FMU's directional
# derivatives only cover the DAE residuals.
_FINITE_DIFFERENCES = {"uses_directional_derivatives": False, "uses_adjoint_derivatives": False}


# Forward/reverse derivatives of an FmuFunction call are computed from its Jacobian
# (jacobian() with coloring, per call). Work-around for a CasADi FmuFunction bug: in a
# multi-direction forward derivative (forward(n), n > 1) a direction without any nonzero seed
# gets the sensitivities of another direction instead of zeros (FmuInternal::eval_ad/eval_fd
# return early without writing osens_), which corrupts NLP Jacobians assembled with several
# directions at a time.
_JACOBIAN_ONLY = {"enable_forward": False, "enable_reverse": False}


def _with_derivative_options(opts):
    der = {k: v for k, v in opts.items() if k not in _JACOBIAN_ONLY and not k.startswith("scheme")}
    return {
        **opts,
        "der_options": der,
        "forward_options": der,
        "reverse_options": der,
        "jacobian_options": der,
    }


def _backward_difference(u, times, history=None):
    """
    Backward finite difference of the timeseries ``u`` on the grid ``times`` (the derivative
    that rtc-tools' collocation uses for a state). At t0: the difference to the last point of
    ``u`` (or else of ``history``) before t0; 0 if there is none.
    """
    times = np.asarray(times, dtype=float)
    values = np.interp(times, u.times, u.values)
    d = np.zeros_like(values)
    d[1:] = np.diff(values) / np.diff(times)
    t_before = u.times[u.times < times[0]]
    if t_before.size:
        t_prev = t_before[-1]
        d[0] = (values[0] - np.interp(t_prev, u.times, u.values)) / (times[0] - t_prev)
    elif history is not None and np.any(history.times < times[0]):
        k = np.flatnonzero(history.times < times[0])[-1]
        d[0] = (values[0] - history.values[k]) / (times[0] - history.times[k])
    return times, d


class _DerivativeInputs(AliasDict):
    """
    Constant inputs, with the FMU input ``der(u)`` of a constant input ``u`` computed on access
    from the timeseries of ``u`` (which a mixin above :py:class:`FMUMixin` in the MRO, such as
    :py:class:`CSVMixin`, adds to this same dictionary after ``FMUMixin.constant_inputs``).
    """

    def __init__(self, relation, derivative_of, times, history):
        super().__init__(relation)
        self.__derivative_of = derivative_of
        self.__times = times
        self.__history = history

    def __contains__(self, key):
        if super().__contains__(key):
            return True
        return key in self.__derivative_of and super().__contains__(self.__derivative_of[key])

    def __getitem__(self, key):
        if super().__contains__(key) or key not in self.__derivative_of:
            return super().__getitem__(key)
        u = self.__derivative_of[key]
        times, d = _backward_difference(
            super().__getitem__(u), self.__times(), self.__history.get(u)
        )
        return Timeseries(times, d)

    def get(self, key, default=None):
        return self[key] if key in self else default


class FMUMixin(OptimizationProblem):
    """
    Adds a model given as an FMI 3.0 Model Exchange FMU with the `fmi-ls-dae
    <https://modelica.github.io/fmi-ls-dae/main/>`_ layered standard (DAE residuals exposed) to
    your optimization problem. It replaces :py:class:`ModelicaMixin`.

    The FMU ``<model_folder>/<fmu_model_name>.fmu`` is imported with CasADi's ``DaeBuilder``
    (CasADi with fmi-ls-dae support: DAE mode is switched on in Configuration Mode). The FMU is
    expected to be exported by OpenModelica with ``--daeMode --simCodeTarget=C+Rust`` (and
    ``--fmiInputDerivatives`` when the model differentiates an input).

    Mapping onto rtc-tools' variables:

    * states: the FMU's continuous states ``x``; derivatives: one symbol ``der(<x>)`` per state,
      fed into the FMU's ``der(x)`` algebraic variables (semi-explicit fmi-ls-dae form);
    * algebraic states: the FMU's other algebraic variables, plus one per FMU output and per
      exposed FMU-internal variable (:py:attr:`fmu_expose_locals`), with the residual
      ``y - F_y(t, x, z, u)``;
    * inputs: an FMU input is a constant input if it has a timeseries
      (:py:meth:`fmu_constant_inputs`), else a control. An FMU input ``der(u)``
      (``--fmiInputDerivatives``) of a constant input ``u`` is a constant input whose
      timeseries is the backward finite difference of ``u`` on the optimization grid;
    * parameters: the FMU's parameters (``causality="parameter"``). FMI 3 lets an importer set
      them only before initialization, so their values are baked into the FMU instances: the
      CasADi functions are rebuilt in :py:meth:`transcribe` whenever :py:meth:`parameters`
      returns other values (e.g. ``theta`` in :py:class:`HomotopyMixin`). Calculated parameters
      are reported by :py:meth:`parameters` but cannot be changed.

    Derivatives: the DAE and initial residuals use the FMU's analytic forward directional
    derivatives (``fmi3GetDirectionalDerivative``); FMU outputs use finite differences (the
    directional derivatives of OpenModelica's FMUs only cover the residuals). There are no
    second derivatives, so IPOPT uses ``hessian_approximation=limited-memory`` by default
    (:py:meth:`solver_options`).

    :cvar fmu_model_name:
        Name of the FMU file (without ``.fmu``). Default is the model name.
    :cvar fmu_expose_locals:
        Expose the FMU's continuous local variables as algebraic states (except names matching
        ModelicaMixin's eliminable-variable pattern ``_name``). ``True``, ``False``, or a list
        of names. Default is ``True``.
    :cvar fmu_fixed_states:
        States whose FMU start value is an initial condition (FMI has no ``fixed`` attribute).
        Default is an empty list.
    :cvar fmu_free_initial_input_derivatives:
        If ``True`` (default), an FMU input ``der(u)`` of a constant input ``u`` is an
        optimization variable fixed by its bounds to the backward finite difference of ``u``
        at every time except t0, where it is free. This matches ModelicaMixin, where such a
        ``u`` is typically aliased to a differentiated state whose derivative at t0 is the free
        ``initial_der`` variable. If ``False``, ``der(u)`` is a constant input, with at t0 the
        backward difference to the data or history before t0, or 0.
    """

    #: Name of the FMU file (without .fmu); default: the model name
    fmu_model_name = None

    #: Expose the FMU's continuous local variables as algebraic states
    fmu_expose_locals = True

    #: States whose FMU start value is an initial condition
    fmu_fixed_states = []

    #: FMU inputs der(u) of constant inputs u: free at t0 (as pymoca's initial derivative)
    fmu_free_initial_input_derivatives = True

    # The collocation residual calls the FMU: it is never affine in the sense of is_affine
    check_collocation_linearity = False
    linear_collocation = False

    # An expanded (SX) FMU call gets its forward derivatives one direction at a time, i.e. one
    # Jacobian evaluation of the FMU call per direction (see _JACOBIAN_ONLY); in MX, one per sweep
    expand_dae_residual = False

    def __init__(self, **kwargs):
        assert "model_folder" in kwargs

        if "model_name" in kwargs:
            model_name = kwargs["model_name"]
        elif hasattr(self, "model_name"):
            model_name = self.model_name
        else:
            model_name = self.__class__.__name__
        fmu_name = self.fmu_model_name or model_name
        fmu_path = os.path.join(kwargs.get("fmu_folder", kwargs["model_folder"]), fmu_name + ".fmu")
        logger.info(f"Loading FMU {fmu_path}.")

        self.__fmu = desc = FMUDescription(fmu_path)
        atexit.register(desc.cleanup)
        self.__fmu_model_name = desc.model_name

        if desc.event_indicators:
            raise NotImplementedError(
                f"FMUMixin: the FMU {fmu_path} has event indicators (e.g. from delay() or "
                "when-clauses). Delays and events are not supported on the FMU route."
            )

        # FMU-internal variables exposed as outputs (algebraic states)
        if self.fmu_expose_locals is True:
            locals_ = [
                v.name
                for v in desc.variables.values()
                if v.causality == "local"
                and v.variability == "continuous"
                and v.type == "Float64"
                and not v.name.startswith("der(")
                and not v.name.startswith("_D_")
                and v.derivative is None
                and not _ELIMINABLE.match(v.name)
            ]
        elif self.fmu_expose_locals:
            locals_ = list(self.fmu_expose_locals)
        else:
            locals_ = []

        # Import with CasADi (once to read the variable sets; the instances used for the
        # optimization are created in __build_functions)
        dae = ca.DaeBuilder(desc.model_name, desc.folder)
        fmu_x, fmu_z, fmu_u = list(dae.x()), list(dae.z()), list(dae.u())
        der_of = {dae.der(s): s for s in fmu_x}
        locals_ = [n for n in locals_ if n not in fmu_x + fmu_z + fmu_u]
        self.__fmu_outputs = list(desc.outputs)
        desc.expose_as_outputs(locals_)
        dae = ca.DaeBuilder(desc.model_name, desc.folder)
        # Outputs that are FMU unknowns (e.g. an algebraic variable) need no output residual
        fmu_y = [n for n in dae.y() if n not in fmu_x + fmu_z + fmu_u]
        for s in fmu_x:
            if dae.der(s) not in fmu_z:
                raise ValueError(
                    f"FMUMixin: der({s}) is not among the FMU's algebraic variables; "
                    "the FMU must be exported in DAE mode (fmi-ls-dae, semi-explicit form)"
                )
        self.__fmu_x, self.__fmu_z, self.__fmu_u, self.__fmu_y = fmu_x, fmu_z, fmu_u, fmu_y
        self.__der_of = der_of

        # Inputs: plain inputs and input derivatives der(u)
        self.__input_derivative_of = {}
        for n in fmu_u:
            m = re.fullmatch(r"der\((.*)\)", n)
            if m and m.group(1) in fmu_u:
                self.__input_derivative_of[n] = m.group(1)
        for n in fmu_u:
            if desc.variables[n].type not in ("Float64", "Float32"):
                raise NotImplementedError(
                    f"FMUMixin: input {n} is of type {desc.variables[n].type}; only Float64 "
                    "inputs are supported on the FMU route"
                )

        constant_input_names = set(self.fmu_constant_inputs(**kwargs))
        for d, u in self.__input_derivative_of.items():
            if u not in constant_input_names:
                raise NotImplementedError(
                    f"FMUMixin: FMU input {d} is the derivative of the control {u}. Only "
                    "derivatives of constant inputs (timeseries) are supported; a control would "
                    "need its collocation finite-difference derivative fed into the FMU."
                )

        # Aliases (FMI <Alias> elements; positive sign only)
        self.__alias_relation = AliasRelation()
        for v in desc.variables.values():
            for a in v.aliases:
                self.__alias_relation.add(v.name, a)

        # MX symbols
        sym = {}
        self.__mx = {}
        self.__mx["time"] = [ca.MX.sym("time")]
        self.__mx["states"] = [sym.setdefault(n, ca.MX.sym(n)) for n in fmu_x]
        self.__mx["derivatives"] = [ca.MX.sym(f"der({n})") for n in fmu_x]
        der_sym = dict(zip(fmu_x, self.__mx["derivatives"], strict=True))
        self.__algebraic_names = [n for n in fmu_z if n not in der_of] + fmu_y
        self.__mx["algebraics"] = [sym.setdefault(n, ca.MX.sym(n)) for n in self.__algebraic_names]
        # Inputs der(u): constant inputs, or controls pinned by their bounds except at t0
        self.__pinned_derivatives = (
            dict(self.__input_derivative_of) if self.fmu_free_initial_input_derivatives else {}
        )
        if not self.fmu_free_initial_input_derivatives:
            constant_input_names |= set(self.__input_derivative_of)
        self.__mx["control_inputs"] = [
            sym.setdefault(n, ca.MX.sym(n)) for n in fmu_u if n not in constant_input_names
        ]
        self.__mx["constant_inputs"] = [
            sym.setdefault(n, ca.MX.sym(n)) for n in fmu_u if n in constant_input_names
        ]
        self.__mx["lookup_tables"] = []
        self.__mx["string_parameters"] = []

        # Parameters: settable ones are symbols; calculated ones are reported values only
        self.__settable_parameters = [
            v.name
            for v in desc.variables.values()
            if v.causality == "parameter"
            and v.variability in ("fixed", "tunable")
            and v.type in ("Float64", "Int32", "Boolean")
        ]
        calculated = [
            v.name
            for v in desc.variables.values()
            if v.causality == "calculatedParameter" and v.type in ("Float64", "Int32", "Boolean")
        ]
        self.__mx["parameters"] = [ca.MX.sym(n) for n in self.__settable_parameters]
        self.__default_parameters = {n: desc.variables[n].start for n in self.__settable_parameters}
        try:
            self.__default_parameters.update(desc.evaluate({}, calculated))
        except Exception as e:
            logger.warning(f"FMUMixin: could not evaluate the calculated parameters: {e}")

        # Residual arguments: z with the der(x) entries fed from the derivative symbols
        self.__args = [
            self.__mx["time"][0],
            ca.vertcat(*self.__mx["states"]),
            ca.vertcat(*[der_sym[der_of[n]] if n in der_of else sym[n] for n in fmu_z]),
            ca.vertcat(*[sym[n] for n in fmu_u]),
        ]
        self.__y_syms = ca.vertcat(*[sym[n] for n in fmu_y])
        self.__y_nominals = None
        self.__baked_parameters = None
        self.__build_functions(self.__default_parameters_settable())

        if logger.getEffectiveLevel() == logging.DEBUG:
            for k in [
                "states",
                "derivatives",
                "algebraics",
                "control_inputs",
                "constant_inputs",
                "parameters",
            ]:
                logger.debug(f"FMUMixin: {k}: {', '.join(v.name() for v in self.__mx[k])}")

        super().__init__(**kwargs)

    # ------------------------------------------------------------------ FMU functions
    def __default_parameters_settable(self):
        return {n: self.__default_parameters[n] for n in self.__settable_parameters}

    def __build_functions(self, parameter_values):
        """(Re)create the FmuFunctions with the given parameter values baked in."""
        desc = self.__fmu
        dae = ca.DaeBuilder(desc.model_name, desc.folder)
        for n, value in parameter_values.items():
            if value != desc.variables[n].start:
                dae.set(n, desc.variables[n].python_type(value))
        # CasADi refuses to evaluate finite differences outside the FMU's min/max. The bounds
        # are imposed by rtc-tools instead (bounds()), which IPOPT may relax slightly.
        for n in self.__fmu_x + self.__fmu_z + self.__fmu_u:
            dae.set_min(n, -np.inf)
            dae.set_max(n, np.inf)
        # Directional derivatives through the FmuFunction's Jacobian (see _JACOBIAN_ONLY)
        opts = {"parallelization": "serial", **_JACOBIAN_ONLY}
        f_alg = dae.create(
            "fmu_dae_residual",
            ["t", "x", "z", "u"],
            ["alg"],
            _with_derivative_options({**opts, **_NO_ADJOINT}),
        )
        residual = [f_alg(*self.__args)]
        if self.__fmu_y:
            time = [v.name for v in desc.variables.values() if v.causality == "independent"]
            scheme = {"t": time, "x": self.__fmu_x, "z": self.__fmu_z, "u": self.__fmu_u}
            scheme["y"] = self.__fmu_y
            f_y = dae.create(
                "fmu_outputs",
                ["t", "x", "z", "u"],
                ["y"],
                _with_derivative_options(
                    {
                        **opts,
                        **_FINITE_DIFFERENCES,
                        "scheme_in": ["t", "x", "z", "u"],
                        "scheme_out": ["y"],
                        "scheme": scheme,
                    }
                ),
            )
            residual.append(self.__y_syms - f_y(*self.__args))
            self.__f_y = f_y
        self.__dae_residual = ca.vertcat(*residual)
        if desc.initial_residuals:
            f_init = dae.create(
                "fmu_initial_residual",
                ["t", "x", "z", "u"],
                ["ialg"],
                _with_derivative_options({**opts, **_NO_ADJOINT}),
            )
            self.__initial_residual = f_init(*self.__args)
        else:
            self.__initial_residual = ca.MX()
        self.__baked_parameters = dict(parameter_values)

    def transcribe(self):
        # Parameter values are baked into the FMU instances: rebuild if they changed
        values = []
        for ensemble_member in range(self.ensemble_size):
            p = self.parameters(ensemble_member)
            for n, default in self.__default_parameters.items():
                if n in self.__settable_parameters or n not in p:
                    continue
                if not _equal(p[n], default):
                    raise ValueError(
                        f"FMUMixin: parameter {n} = {p[n]} differs from its value {default} in "
                        "the FMU, where it is a calculated or structural parameter that cannot be "
                        "set. Export the FMU with the desired value instead."
                    )
            v = {}
            for n in self.__settable_parameters:
                value = p.get(n, self.__default_parameters[n])
                if isinstance(value, ca.MX):
                    if not value.is_constant():
                        raise ValueError(f"FMUMixin: parameter {n} must be numeric, got {value}")
                    value = float(value)
                v[n] = value
            values.append(v)
        if any(any(not _equal(v[n], values[0][n]) for n in v) for v in values[1:]):
            raise NotImplementedError(
                "FMUMixin: ensemble members with different FMU parameter values"
            )
        if any(not _equal(values[0][n], self.__baked_parameters[n]) for n in values[0]):
            changed = {
                n: values[0][n]
                for n in values[0]
                if not _equal(values[0][n], self.__baked_parameters[n])
            }
            logger.info(f"FMUMixin: new FMU instances for parameter values {changed}")
            self.__build_functions(values[0])
            self.clear_transcription_cache()
        if self.__y_nominals is None:
            self.__y_nominals = self.__estimate_output_nominals()
        return super().transcribe()

    def __estimate_output_nominals(self):
        """
        Nominals of the FMU outputs (and exposed FMU-internal variables) that have no nominal
        attribute in the FMU: ``max(|y|, sum_j |dy/dv_j| nominal_j)`` at a reference point
        (history, seed or constant input at t0, else the midpoint of finite bounds, else the
        FMU start value), if that exceeds 1. ModelicaMixin gets these from pymoca's alias
        elimination (e.g. ``V_storage = storage.V`` has the nominal of ``storage.V``).
        """
        missing = [n for n in self.__fmu_y if self.__fmu.variables[n].nominal == 1.0]
        if not missing:
            return {}
        t0 = self.initial_time
        history, seed = self.history(0), self.seed(0)
        constant_inputs, bounds = self.constant_inputs(0), self.bounds()

        def at_t0(ts):
            try:
                return float(np.interp(t0, ts.times, ts.values))
            except Exception:
                return None

        def reference(n):
            for d in (history, seed, constant_inputs):
                if n in d:
                    v = at_t0(d[n])
                    if v is not None and np.isfinite(v):
                        return v
            m, M = bounds.get(n, (-np.inf, np.inf))
            if not isinstance(m, Timeseries) and not isinstance(M, Timeseries):
                if np.isfinite(m) and np.isfinite(M):
                    return 0.5 * (m + M)
            start = self.__fmu.variables[n].start
            return float(start) if start is not None else 0.0

        x = np.array([reference(n) for n in self.__fmu_x])
        z = np.array([0.0 if n in self.__der_of else reference(n) for n in self.__fmu_z])
        u = np.array([reference(n) for n in self.__fmu_u])
        nominals = np.concatenate(
            [
                [self.variable_nominal(n) for n in self.__fmu_x],
                [1.0 if n in self.__der_of else self.variable_nominal(n) for n in self.__fmu_z],
                [self.variable_nominal(n) for n in self.__fmu_u],
            ]
        )
        try:
            y = np.array(self.__f_y(0.0, x, z, u)).ravel()
            jac = self.__f_y.jacobian()(0.0, x, z, u, y)
            J = np.hstack([np.array(ca.densify(j)) for j in jac[1:]])
        except Exception as e:
            logger.warning(f"FMUMixin: could not estimate nominals of the FMU outputs: {e}")
            return {}
        estimate = np.maximum(np.abs(y), np.abs(J) @ nominals)
        result = {}
        for i, n in enumerate(self.__fmu_y):
            if n in missing and np.isfinite(estimate[i]) and estimate[i] > 1.0:
                result[n] = float(estimate[i])
                logger.debug(f"FMUMixin: estimated nominal of {n}: {result[n]}")
        return result

    # ------------------------------------------------------------------ rtc-tools interface
    def fmu_constant_inputs(self, **kwargs) -> list[str]:
        """
        Names of the FMU inputs that are constant inputs (given by data) rather than controls.

        ModelicaMixin uses the inputs' ``fixed`` attribute, which FMI does not have. The
        default implementation takes the inputs that have a column in the CSV timeseries import
        (``<input_folder>/<timeseries_import_basename>.csv``, or that of the first ensemble
        member). Override it for other input formats.
        """
        input_folder = kwargs.get("input_folder")
        basename = getattr(self, "timeseries_import_basename", "timeseries_import")
        candidates = [os.path.join(input_folder or "", basename + ".csv")]
        if input_folder and os.path.isdir(input_folder):
            candidates += [
                os.path.join(input_folder, d, basename + ".csv")
                for d in sorted(os.listdir(input_folder))
            ]
        delimiter = getattr(self, "csv_delimiter", ",")
        for path in candidates:
            if os.path.isfile(path):
                with open(path, newline="") as f:
                    header = next(csv.reader(f, delimiter=delimiter), [])
                columns = {c.strip() for c in header[1:]}
                return [
                    n
                    for n in self.__fmu_u
                    if n in columns or any(a in columns for a in self.__fmu.variables[n].aliases)
                ]
        logger.warning(
            "FMUMixin: no CSV timeseries import found; all FMU inputs are controls. "
            "Override fmu_constant_inputs() to declare constant inputs."
        )
        return []

    @property
    def dae_residual(self):
        return self.__dae_residual

    @property
    def initial_residual(self):
        return self.__initial_residual

    @property
    def dae_variables(self):
        return self.__mx

    @property
    @cached
    def alias_relation(self):
        return self.__alias_relation

    @property
    @cached
    def output_variables(self):
        output_variables = [ca.MX.sym(n) for n in self.__fmu_outputs]
        output_variables.extend(
            v for v in self.__mx["control_inputs"] if v.name() not in self.__pinned_derivatives
        )
        return output_variables

    def map_options(self):
        # One FMU instance serves all time steps; no parallel evaluation of the FMU
        options = super().map_options()
        options["mode"] = "serial"
        return options

    def solver_options(self):
        options = super().solver_options()
        if options.get("solver") == "ipopt":
            # The FMU provides no second derivatives
            options["ipopt"].setdefault("hessian_approximation", "limited-memory")
        # Each forward sweep through an FMU call evaluates that call's full Jacobian (see
        # _JACOBIAN_ONLY): many directions per sweep in the DAE functions (max_num_dir)
        options["optimized_num_dir"] = max(options.get("optimized_num_dir", 0), 64)
        return options

    @cached
    def parameters(self, ensemble_member):
        parameters = super().parameters(ensemble_member)
        parameters.update(self.__default_parameters)
        return parameters

    @cached
    def constant_inputs(self, ensemble_member):
        constant_inputs = super().constant_inputs(ensemble_member)
        if not self.__input_derivative_of or self.__pinned_derivatives:
            return constant_inputs
        lazy = _DerivativeInputs(
            self.alias_relation,
            self.__input_derivative_of,
            self.times,
            self.history(ensemble_member),
        )
        lazy.update(constant_inputs)
        return lazy

    @cached
    def history(self, ensemble_member):
        history = super().history(ensemble_member)
        initial_time = np.array([self.initial_time])
        for n in self.fmu_fixed_states:
            start = self.__fmu.variables[n].start
            if start is None:
                raise ValueError(f"FMUMixin: fixed state {n} has no start value in the FMU")
            history[n] = Timeseries(initial_time, float(start))
        return history

    def __model_variables(self):
        return self.__fmu_x + self.__algebraic_names + self.__fmu_u

    @cached
    @ensemble_bounds_check
    def bounds(self, ensemble_member: int | None = None):
        bounds = (
            super().bounds(ensemble_member) if self.ensemble_specific_bounds else super().bounds()
        )
        for n in self.__model_variables():
            v = self.__fmu.variables[n]
            m, M = bounds.get(n, (-np.inf, np.inf))
            bounds[n] = (max(m, v.min), min(M, v.max))
        for d in self.__pinned_derivatives:
            times, values = self.__pinned_derivative(d, ensemble_member or 0)
            lower, upper = values.copy(), values.copy()
            lower[0], upper[0] = -np.inf, np.inf
            bounds[d] = (Timeseries(times, lower), Timeseries(times, upper))
        return bounds

    def __pinned_derivative(self, d, ensemble_member):
        u = self.__pinned_derivatives[d]
        try:
            ts = self.constant_inputs(ensemble_member)[u]
        except KeyError:
            raise Exception(f"FMUMixin: no values found for constant input {u} (for {d})")
        times, values = _backward_difference(ts, self.times())
        values[0] = 0.0
        return times, values

    @cached
    def seed(self, ensemble_member):
        seed = super().seed(ensemble_member)
        for d in self.__pinned_derivatives:
            seed[d] = Timeseries(*self.__pinned_derivative(d, ensemble_member))
        for n in self.__fmu_x + self.__algebraic_names:
            if n in self.fmu_fixed_states:
                continue
            start = self.__fmu.variables[n].start
            if start is not None and start != 0.0 and np.isfinite(start):
                times = self.times(n)
                seed[n] = Timeseries(times, np.full_like(times, float(start)))
        return seed

    def variable_is_discrete(self, variable):
        return False

    @property
    @cached
    def __nominals(self):
        nominals = AliasDict(self.alias_relation)
        for n in self.__model_variables():
            nominal = abs(self.__fmu.variables[n].nominal)
            if np.isfinite(nominal) and nominal not in (0.0, 1.0):
                nominals[n] = nominal
        return nominals

    def variable_nominal(self, variable):
        if self.__y_nominals:
            canonical, _ = self.alias_relation.canonical_signed(variable)
            if canonical in self.__y_nominals:
                return self.__y_nominals[canonical]
        try:
            return self.__nominals[variable]
        except KeyError:
            return super().variable_nominal(variable)


def _equal(a, b):
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return a == b
