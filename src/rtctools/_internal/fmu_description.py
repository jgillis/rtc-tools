"""
Reading an FMI 3.0 Model Exchange FMU with the fmi-ls-dae layered standard.

CasADi's ``DaeBuilder`` imports the FMU, but it does not expose everything rtc-tools needs:
FMI ``<Alias>`` elements, the causality/variability/initial attributes of every variable, and
the values of calculated parameters. This module reads them from ``modelDescription.xml`` and
the fmi-ls-dae manifest, and evaluates parameter values through the FMI C API.
"""

import ctypes
import glob
import os
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass, field

import defusedxml.ElementTree as DefusedElementTree

LS_DAE_MANIFEST = "extra/org.fmi-standard.fmi-ls-dae/fmi-ls-manifest.xml"
INITIAL_RESIDUALS_ANNOTATION = "org.openmodelica.fmi-ls-dae.initial-residuals"

_TYPES = {"Float64": float, "Float32": float, "Boolean": bool}


@dataclass
class FMUVariable:
    name: str
    value_reference: int
    type: str
    causality: str = "local"
    variability: str = "continuous"
    initial: str | None = None
    start: float | bool | int | None = None
    min: float = float("-inf")
    max: float = float("inf")
    nominal: float = 1.0
    description: str = ""
    aliases: list[str] = field(default_factory=list)
    derivative: int | None = None

    @property
    def python_type(self):
        return _TYPES.get(self.type, int)

    @property
    def is_parameter(self):
        return self.causality in ("parameter", "calculatedParameter")


def _parse_value(type_, text):
    if text is None:
        return None
    if type_ == "Boolean":
        return text.strip().lower() in ("true", "1")
    if type_ in ("Float64", "Float32"):
        return float(text)
    if type_.startswith("Int") or type_.startswith("UInt") or type_ == "Enumeration":
        return int(text)
    return text


class FMUDescription:
    """
    The variables of an FMI 3.0 FMU (``modelDescription.xml``) and its fmi-ls-dae manifest.

    The FMU is extracted to a temporary folder (:py:attr:`folder`), which CasADi's
    ``DaeBuilder`` can load directly.
    """

    def __init__(self, fmu_path: str):
        self.fmu_path = fmu_path
        if os.path.isdir(fmu_path):
            src = fmu_path
            self.folder = tempfile.mkdtemp(prefix="rtctools_fmu_")
            shutil.copytree(src, self.folder, dirs_exist_ok=True)
        else:
            self.folder = tempfile.mkdtemp(prefix="rtctools_fmu_")
            with zipfile.ZipFile(fmu_path) as z:
                z.extractall(self.folder)

        self.__md_path = os.path.join(self.folder, "modelDescription.xml")
        md_root = DefusedElementTree.parse(self.__md_path).getroot()
        self.fmi_version = md_root.get("fmiVersion", "")
        if not self.fmi_version.startswith("3."):
            raise ValueError(f"FMU {fmu_path}: FMI {self.fmi_version} is not supported (needs 3.0)")
        self.model_name = md_root.get("modelName")
        self.instantiation_token = md_root.get("instantiationToken")
        me = md_root.find("ModelExchange")
        if me is None:
            raise ValueError(f"FMU {fmu_path} is not a Model Exchange FMU")
        self.model_identifier = me.get("modelIdentifier")
        self.provides_directional_derivatives = (
            me.get("providesDirectionalDerivatives", "false") == "true"
        )

        self.variables: dict[str, FMUVariable] = {}
        self.by_vr: dict[int, FMUVariable] = {}
        self.alias_of: dict[str, str] = {}
        for v in md_root.find("ModelVariables"):
            type_ = v.tag
            var = FMUVariable(
                name=v.get("name"),
                value_reference=int(v.get("valueReference")),
                type=type_,
                causality=v.get("causality", "local"),
                variability=v.get(
                    "variability",
                    "continuous" if type_ in ("Float64", "Float32") else "discrete",
                ),
                initial=v.get("initial"),
                start=_parse_value(type_, v.get("start")),
                description=v.get("description", ""),
            )
            for attr in ("min", "max", "nominal"):
                if v.get(attr) is not None and type_ in ("Float64", "Float32"):
                    setattr(var, attr, float(v.get(attr)))
            if v.get("derivative") is not None:
                var.derivative = int(v.get("derivative"))
            for a in v.findall("Alias"):
                var.aliases.append(a.get("name"))
                self.alias_of[a.get("name")] = var.name
            self.variables[var.name] = var
            self.by_vr[var.value_reference] = var

        ms = md_root.find("ModelStructure")
        self.event_indicators = (
            [] if ms is None else [int(e.get("valueReference")) for e in ms.findall("EventIndicator")]
        )
        self.outputs = (
            [] if ms is None else [self.by_vr[int(e.get("valueReference"))].name
                                   for e in ms.findall("Output")]
        )

        # fmi-ls-dae manifest
        mf_path = os.path.join(self.folder, LS_DAE_MANIFEST)
        if not os.path.exists(mf_path):
            raise ValueError(
                f"FMU {fmu_path} has no fmi-ls-dae manifest ({LS_DAE_MANIFEST}); "
                "export it with OpenModelica --daeMode"
            )
        mf_root = DefusedElementTree.parse(mf_path).getroot()
        self.dae_residuals = []
        self.initial_residuals = []
        for r in mf_root.iter("Residual"):
            name = self.by_vr[int(r.get("valueReference"))].name
            if name.startswith("_D_INITres"):
                self.initial_residuals.append(name)
            else:
                self.dae_residuals.append(name)
        enable = mf_root.find(".//EnableDAEParameter")
        self.dae_mode_parameter = (
            self.by_vr[int(enable.get("valueReference"))].name
            if enable is not None and enable.get("valueReference") is not None
            else ("_D_daeMode" if "_D_daeMode" in self.variables else None)
        )

    def cleanup(self):
        shutil.rmtree(self.folder, ignore_errors=True)

    def expose_as_outputs(self, names: list[str]):
        """
        Mark FMU-internal (local) variables as outputs in the extracted modelDescription.xml, with
        ``<Output>`` elements without ``dependencies`` (i.e. depending on all knowns).

        CasADi only knows the dependencies of outputs, residuals and state derivatives; for a
        local variable it assumes there are none, so that its Jacobian would be zero. The FMU
        binary does not read modelDescription.xml, so this changes nothing in the evaluation.
        """
        if not names:
            return
        text = open(self.__md_path, encoding="utf-8").read()
        new_outputs = []
        for n in names:
            var = self.variables[n]
            pat = re.compile(r'(<\w+\s+name="' + re.escape(n) + r'"\s+valueReference="'
                             + str(var.value_reference) + r'")')
            text, count = pat.subn(r'\1 causality="output"', text, count=1)
            if count != 1:
                raise RuntimeError(f"Cannot expose FMU variable {n} as an output")
            # An explicit causality="local" would now be duplicated
            text = re.sub(r'(name="' + re.escape(n) + r'"[^>]*?causality="output")'
                          r'([^>]*?)\s+causality="local"', r"\1\2", text, count=1)
            var.causality = "output"
            new_outputs.append(f'    <Output valueReference="{var.value_reference}"/>\n')
        open(self.__md_path, "w", encoding="utf-8").write(_add_outputs(text, new_outputs))
        # The manifest's ModelStructure replaces the Output elements of modelDescription.xml
        mf_path = os.path.join(self.folder, LS_DAE_MANIFEST)
        mf_text = open(mf_path, encoding="utf-8").read()
        if "<ModelStructure" in mf_text:
            open(mf_path, "w", encoding="utf-8").write(_add_outputs(mf_text, new_outputs))
        self.outputs = names_unique(self.outputs + list(names))

    def evaluate(self, start_values: dict, names: list[str]) -> dict:
        """
        Instantiate the FMU through the FMI 3 C API, switch DAE mode on, set the given start
        values (parameters), run Initialization Mode and read the variables in ``names`` (e.g.
        calculated parameters, which have no start value in modelDescription.xml).
        """
        libs = glob.glob(os.path.join(self.folder, "binaries", "*", self.model_identifier + ".*"))
        libs = [p for p in libs if p.endswith((".so", ".dylib", ".dll"))]
        if not libs:
            raise RuntimeError(f"No FMU binary for {self.model_identifier} in {self.folder}")
        lib = ctypes.CDLL(libs[0])
        CB = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p)
        log_messages = []
        cb = CB(lambda env, status, cat, msg: log_messages.append(msg.decode(errors="replace")))
        lib.fmi3InstantiateModelExchange.restype = ctypes.c_void_p
        resources = (os.path.join(self.folder, "resources") + os.sep).encode()
        inst = lib.fmi3InstantiateModelExchange(
            b"rtctools_parameters", self.instantiation_token.encode(), resources, False, False,
            None, cb,
        )
        if not inst:
            raise RuntimeError("fmi3InstantiateModelExchange failed")
        inst = ctypes.c_void_p(inst)

        ctype = {"Float64": ctypes.c_double, "Int32": ctypes.c_int32, "Boolean": ctypes.c_bool}

        def call(fname, *args):
            status = getattr(lib, fname)(inst, *args)
            if status > 1:  # fmi3Error, fmi3Fatal
                raise RuntimeError(f"{fname} failed: {'; '.join(log_messages[-3:])}")

        def set_value(var, value, setter="Set"):
            ct = ctype.get(var.type)
            if ct is None:
                raise NotImplementedError(f"FMU variable {var.name} of type {var.type}")
            call(f"fmi3{setter}{var.type}", (ctypes.c_uint32 * 1)(var.value_reference),
                 ctypes.c_size_t(1), (ct * 1)(value), ctypes.c_size_t(1))

        try:
            if self.dae_mode_parameter is not None:
                call("fmi3EnterConfigurationMode")
                set_value(self.variables[self.dae_mode_parameter], True)
                call("fmi3ExitConfigurationMode")
            for n, value in start_values.items():
                set_value(self.variables[n], value)
            call("fmi3EnterInitializationMode", ctypes.c_bool(False), ctypes.c_double(0.0),
                 ctypes.c_double(0.0), ctypes.c_bool(False), ctypes.c_double(0.0))
            # Read in Initialization Mode, variable by variable: a Float64 get may evaluate the
            # model equations, which can fail at the default (start) values.
            result = {}
            for n in names:
                var = self.variables[n]
                ct = ctype.get(var.type)
                if ct is None:
                    continue
                buf = (ct * 1)()
                status = getattr(lib, f"fmi3Get{var.type}")(
                    inst, (ctypes.c_uint32 * 1)(var.value_reference), ctypes.c_size_t(1), buf,
                    ctypes.c_size_t(1),
                )
                if status <= 1:
                    result[n] = var.python_type(buf[0])
            return result
        finally:
            lib.fmi3FreeInstance(inst)


def names_unique(names):
    seen = set()
    return [n for n in names if not (n in seen or seen.add(n))]


def _add_outputs(text, outputs):
    """Insert <Output> elements at the start of the (first) ModelStructure element."""
    if "<ModelStructure>" in text:
        return text.replace("<ModelStructure>", "<ModelStructure>\n" + "".join(outputs), 1)
    return text.replace(
        "<ModelStructure/>", "<ModelStructure>\n" + "".join(outputs) + "  </ModelStructure>", 1
    )
