# RTC-Tools with the Base Modelica and FMU routes

Work in progress, shared to try out. This image puts together changes to several projects that let
RTC-Tools use a model compiled by OpenModelica, instead of only pymoca's own Modelica front end:

- **Base Modelica route:** OpenModelica flattens the model to Base Modelica, and pymoca reads that.
- **FMU route:** OpenModelica exports an FMI 3.0 FMU with the fmi-ls-dae layered standard, and the
  new `FMUMixin` uses it as the model through CasADi.

The changes are on their way upstream in separate pull requests. Until then, each project has a
branch `rtc-basemodelica` on a fork that combines them:

| Project | Branch | Contents |
|---|---|---|
| [rtc-tools](https://github.com/jgillis/rtc-tools/tree/rtc-basemodelica) | `rtc-basemodelica` | Example models conforming to the Modelica specification, `FMUMixin` and FMU drivers, OpenModelica conformance test, this folder |
| [pymoca](https://github.com/jgillis/pymoca/tree/rtc-basemodelica) | `rtc-basemodelica` | Array comprehensions, inverse trigonometric functions, `A[i, i]` in loops, `basemodelica` compiler option |
| [OpenModelica](https://github.com/jgillis/OpenModelica/tree/rtc-basemodelica) | `rtc-basemodelica` | fmi-ls-dae export with initial residuals, analytic derivatives and `--fmiInputDerivatives` |
| [CasADi](https://github.com/jgillis/casadi/tree/rtc-basemodelica) | `rtc-basemodelica` | Enabling fmi-ls-dae in Configuration Mode, initial residuals |
| [rtc-tools-channel-flow](https://github.com/jgillis/rtc-tools-channel-flow/tree/rtc-basemodelica) | `rtc-basemodelica` | Valid Modelica for the `H.min` bound in `PartialHomotopic` |
| [fmi-ls-dae](https://github.com/jgillis/fmi-ls-dae/tree/rtc-basemodelica) | `rtc-basemodelica` | Specification text for initialization in DAE mode (documentation only, not used by the build) |

## Build

```bash
docker build --pull -t rtc-basemodelica "https://github.com/jgillis/rtc-tools.git#rtc-basemodelica:docker/rtc-basemodelica"
```

OpenModelica and CasADi are compiled from source, so a first build takes a while and needs about
25 GB of free disk. Use an x86-64 machine: other platforms have not been tried.

## Use

```bash
docker run --rm -it rtc-basemodelica
```

You start in `/src/rtc-tools`. The Python packages in `/src` are editable installs, and the examples
already have their exports next to their models (`*.basemodelica.mo`, `*.fmu`).

Run an example on each route:

```bash
cd /src/rtc-tools/examples/basic/src
python example.py                                                    # pymoca, as usual
PYTHONPATH=/opt/rtc-basemodelica:$PYTHONPATH python example.py       # Base Modelica
python example_fmu.py                                                # FMU
```

Run the tests:

```bash
/opt/rtc-basemodelica/run_tests.sh                # everything
/opt/rtc-basemodelica/run_tests.sh unit omc       # some parts: unit examples examples-bm examples-fmu omc
```

To export again after changing a model, run `/opt/rtc-basemodelica/export_basemodelica.sh` or
`/opt/rtc-basemodelica/export_fmus.sh`.

## Known limitations

- The Base Modelica route is switched on from outside with `sitecustomize.py`: `ModelicaMixin` has
  no option for it yet.
- On the FMU route, `FMUMixin` uses a limited-memory Hessian, which makes large problems much slower
  than pymoca with its exact Hessian (`channel_wave_damping/example_optimization_fmu.py` takes about
  6 minutes, against about 20 seconds).
- The FMU route has drivers for 6 scripts of 5 examples. Examples with delays, events or Boolean
  variables do not run on it yet.
- `fmi3GetAdjointDerivative` of the exported FMU gives wrong results on `channel_pulse`, so
  `FMUMixin` does not rely on it.
- The drivers do not compare the results of the routes: passing here means that the script ran to
  the end.
- `tests/modelica/test_omc_conformance.py` still records known messages. `cascading_channels` and
  `pumped_hydropower_system` do not translate in OpenModelica (under-determined), and
  `channel_wave_damping` translates with warnings about its initial conditions.
