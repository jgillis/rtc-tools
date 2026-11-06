#!/bin/bash
# Export the example models as Base Modelica to examples/<example>/model/<class>.basemodelica.mo,
# which pymoca reads with the basemodelica compiler option (run_tests.sh examples-bm).
set -euo pipefail

RTC=${RTC_TOOLS:-/src/rtc-tools}
OMC=${OMC:-omc}
MSL=${MSL:-$HOME/.openmodelica/libraries/Modelica 4.0.0+maint.om}
CHANNEL_FLOW=${CHANNEL_FLOW:-$(python -c 'import os, rtctools_channel_flow as m; print(os.path.dirname(m.__file__))')/modelica}

export_basemodelica() {
  local example=$1 class=$2
  local model_dir=$RTC/examples/$example/model
  local work file
  work=$(mktemp -d)
  {
    echo "loadFile(\"$MSL/package.mo\"); getErrorString();"
    echo "loadFile(\"$CHANNEL_FLOW/Deltares/package.mo\"); getErrorString();"
    for file in "$model_dir"/*.mo; do
      case $file in
        *.basemodelica.mo) ;;
        *) echo "loadFile(\"$file\"); getErrorString();" ;;
      esac
    done
    echo "writeFile(\"$model_dir/$class.basemodelica.mo\", instantiateModel($class)); getErrorString();"
  } > "$work/export.mos"
  (cd "$work" && "$OMC" --BaseModelica --evaluateStructuralParameters=strictlyNecessary \
      --baseModelicaFormat=scalarized --baseModelicaOptions=scalarize,moveBindings export.mos > omc.log 2>&1) || true
  if [ -s "$model_dir/$class.basemodelica.mo" ]; then
    echo "$example/model/$class.basemodelica.mo"
  else
    echo "FAILED: $example/model/$class.basemodelica.mo" >&2
    tail -30 "$work/omc.log" >&2
    rm -rf "$work"
    rm -f "$model_dir/$class.basemodelica.mo"
    return 1
  fi
  rm -rf "$work"
}

export_basemodelica basic Example
export_basemodelica cascading_channels Example
export_basemodelica channel_pulse Example
export_basemodelica channel_wave_damping ExampleLocalControl
export_basemodelica channel_wave_damping ExampleOptimization
export_basemodelica ensemble Example
export_basemodelica goal_programming Example
export_basemodelica integrator_delay Example
export_basemodelica lookup_table Example
export_basemodelica mixed_integer Example
export_basemodelica pumped_hydropower_system PumpedStoragePlant
export_basemodelica simulation Example
export_basemodelica step_response Example
