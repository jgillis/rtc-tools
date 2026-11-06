"""Compile with pymoca's basemodelica option (run_tests.sh examples-bm).

RTC-Tools has no switch for this route yet: this sets the option on both mixins from outside.
"""

import rtctools.optimization.modelica_mixin as _modelica_mixin
import rtctools.simulation.simulation_problem as _simulation_problem

for _cls in (_modelica_mixin.ModelicaMixin, _simulation_problem.SimulationProblem):

    def _compiler_options(self, _original=_cls.compiler_options):
        options = _original(self)
        options["basemodelica"] = True
        return options

    _cls.compiler_options = _compiler_options
