#!/bin/bash
# Export the FMUs of the FMU route drivers (examples/*/src/*_fmu.py) to examples/<example>/model/.
#
# FMI 3.0 Model Exchange with the fmi-ls-dae layered standard. OpenModelica writes only literal
# min/max/nominal, so where those depend on a parameter the exported class is a wrapper that
# restates them. stateSelect makes the states H and Q, as with pymoca, instead of internal
# variables such as Channel._wetted_perimeter. Parameters that OpenModelica evaluates as structural
# (use_convective_acceleration, use_upwind) need one FMU per value.
set -euo pipefail

RTC=${RTC_TOOLS:-/src/rtc-tools}
OMC=${OMC:-omc}
MSL=${MSL:-$HOME/.openmodelica/libraries/Modelica 4.0.0+maint.om}
CHANNEL_FLOW=${CHANNEL_FLOW:-$(python -c 'import os, rtctools_channel_flow as m; print(os.path.dirname(m.__file__))')/modelica}
OPTIONS="--daeMode --simCodeTarget=C+Rust --fmiInputDerivatives"

# export_fmu <example> <output name> <class> <wrapper or ""> <model file>...
export_fmu() {
  local example=$1 output=$2 class=$3 wrapper=$4
  shift 4
  local model_dir=$RTC/examples/$example/model
  local work
  work=$(mktemp -d)
  {
    echo "loadFile(\"$MSL/package.mo\"); getErrorString();"
    echo "loadFile(\"$CHANNEL_FLOW/Deltares/package.mo\"); getErrorString();"
    for file in "$@"; do
      echo "loadFile(\"$model_dir/$file\"); getErrorString();"
    done
    if [ -n "$wrapper" ]; then
      echo "$wrapper" > "$work/wrapper.mo"
      echo "loadFile(\"$work/wrapper.mo\"); getErrorString();"
    fi
    echo "setCommandLineOptions(\"$OPTIONS\"); getErrorString();"
    echo "buildModelFMU($class, version=\"3.0\", fmuType=\"me\", fileNamePrefix=\"$class\"); getErrorString();"
  } > "$work/export.mos"
  if (cd "$work" && "$OMC" export.mos > omc.log 2>&1) && [ -f "$work/$class.fmu" ]; then
    cp "$work/$class.fmu" "$model_dir/$output.fmu"
    echo "$example/model/$output.fmu"
  else
    echo "FAILED: $example/model/$output.fmu" >&2
    tail -30 "$work/omc.log" >&2
    rm -rf "$work"
    return 1
  fi
  rm -rf "$work"
}

export_fmu basic Example Example "" Example.mo
export_fmu single_reservoir SingleReservoir SingleReservoir "" single_reservoir.mo

export_fmu step_response Example Example_fmu \
  'model Example_fmu extends Example(Channel(H(each min = 0.2), Q(each nominal = 100))); end Example_fmu;' \
  Example.mo
export_fmu step_response Example_convective Example_fmu_convective \
  'model Example_fmu_convective extends Example(Channel(use_convective_acceleration = true, H(each min = 0.2), Q(each nominal = 100))); end Example_fmu_convective;' \
  Example.mo

export_fmu channel_pulse Example Example_fmu \
  'model Example_fmu extends Example(Channel(H(each min = -5.0, each stateSelect = StateSelect.prefer), Q(each nominal = 100))); end Example_fmu;' \
  Example.mo
export_fmu channel_pulse Example_convective Example_fmu_convective \
  'model Example_fmu_convective extends Example(Channel(use_convective_acceleration = true, H(each min = -5.0, each stateSelect = StateSelect.prefer), Q(each nominal = 100))); end Example_fmu_convective;' \
  Example.mo
export_fmu channel_pulse Example_upwind Example_fmu_upwind \
  'model Example_fmu_upwind extends Example(Channel(use_convective_acceleration = true, use_upwind = true, H(each min = -5.0, each stateSelect = StateSelect.prefer), Q(each nominal = 100))); end Example_fmu_upwind;' \
  Example.mo

reaches='upstream(H(each min = 15, each stateSelect = StateSelect.prefer), Q(each nominal = 100)), middle(H(each min = 10, each stateSelect = StateSelect.prefer), Q(each nominal = 100)), downstream(H(each min = 5, each stateSelect = StateSelect.prefer), Q(each nominal = 100))'
export_fmu channel_wave_damping ExampleOptimization ExampleOptimization_fmu \
  "model ExampleOptimization_fmu extends ExampleOptimization($reaches); end ExampleOptimization_fmu;" \
  Example.mo ExampleOptimization.mo PIDController.mo
export_fmu channel_wave_damping ExampleLocalControl ExampleLocalControl_fmu \
  "model ExampleLocalControl_fmu extends ExampleLocalControl($reaches); end ExampleLocalControl_fmu;" \
  Example.mo ExampleLocalControl.mo PIDController.mo
