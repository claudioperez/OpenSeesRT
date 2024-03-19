"""
This module implements the OpenSeesPy interface.
Imports can be performed exactly as one would
from openseespy, for example:

>>> import opensees.openseespy as ops

>>> from opensees.openseespy import node, model

>>> from opensees.openseespy import *


The fiber section API is one of the only unsupported
features. This may be implemented in the
near future, but even then, using it would be strongly
discouraged.
"""
import json
from functools import partial

from .tcl import Interpreter

# A list of symbol names that are importable
# from this module. All of these are dynamically
# resolved by the function __getattr__ below.
__all__ = [
    "uniaxialMaterial",
    "testUniaxialMaterial",
    "setStrain",
    "getStrain",
    "getStress",
    "getTangent",
    "getDampTangent",
    "wipe",
    "model",
    "node",
    "fix",
    "element",
    "timeSeries",
    "pattern",
    "load",
    "system",
    "numberer",
    "constraints",
    "integrator",
    "algorithm",
    "analysis",
    "analyze",
    "test",
    "section",
    "fiber",
    "patch",
    "layer",
    "geomTransf",
    "beamIntegration",
    "loadConst",
    "eleLoad",
    "reactions",
    "nodeReaction",
    "eigen",
    "modalProperties",
    "responseSpectrumAnalysis",
    "nDMaterial",
    "block2D",
    "block3D",
    "rayleigh",
    "wipeAnalysis",
    "setTime",
    "remove",
    "mass",
    "equalDOF",
    "nodeEigenvector",
    "getTime",
    "setCreep",
    "eleResponse",
    "sp",
    "fixX",
    "fixY",
    "fixZ",
    "reset",
    "initialize",
    "getLoadFactor",
    "build",
    "printModel",
    "printA",
    "printB",
    "printGID",
    "testNorm",
    "testNorms",
    "testIter",
    "recorder",
    "database",
    "save",
    "restore",
    "eleForce",
    "eleDynamicalForce",
    "nodeUnbalance",
    "nodeDisp",
    "setNodeDisp",
    "nodeVel",
    "setNodeVel",
    "nodeAccel",
    "setNodeAccel",
    "nodeResponse",
    "nodeCoord",
    "setNodeCoord",
    "getPatterns",
    "getFixedNodes",
    "getFixedDOFs",
    "getConstrainedNodes",
    "getConstrainedDOFs",
    "getRetainedNodes",
    "getRetainedDOFs",
    "updateElementDomain",
    "getNDM",
    "getNDF",
    "eleNodes",
    "eleType",
    "nodeDOFs",
    "nodeMass",
    "nodePressure",
    "setNodePressure",
    "nodeBounds",
    "start",
    "stop",
    "modalDamping",
    "modalDampingQ",
    "setElementRayleighDampingFactors",
    "region",
    "setPrecision",
    "searchPeerNGA",
    "domainChange",
    "record",
    "metaData",
    "defaultUnits",
    "stripXML",
    "convertBinaryToText",
    "convertTextToBinary",
    "getEleTags",
    "getCrdTransfTags",
    "getNodeTags",
    "getParamTags",
    "getParamValue",
    "sectionForce",
    "sectionDeformation",
    "sectionStiffness",
    "sectionFlexibility",
    "sectionLocation",
    "sectionWeight",
    "sectionTag",
    "sectionDisplacement",
    "cbdiDisplacement",
    "basicDeformation",
    "basicForce",
    "basicStiffness",
    "InitialStateAnalysis",
    "totalCPU",
    "solveCPU",
    "accelCPU",
    "numFact",
    "numIter",
    "systemSize",
    "version",
    "setMaxOpenFiles",
    "limitCurve",
    "imposedMotion",
    "imposedSupportMotion",
    "groundMotion",
    "equalDOF_Mixed",
    "rigidLink",
    "rigidDiaphragm",
    "ShallowFoundationGen",
    "setElementRayleighFactors",
    "mesh",
    "remesh",
    "parameter",
    "addToParameter",
    "updateParameter",
    "setParameter",
    "getPID",
    "getNP",
    "barrier",
    "send",
    "recv",
    "Bcast",
    "frictionModel",
    "computeGradients",
    "sensitivityAlgorithm",
    "sensNodeDisp",
    "sensNodeVel",
    "sensNodeAccel",
    "sensLambda",
    "sensSectionForce",
    "sensNodePressure",
    "getNumElements",
    "getEleClassTags",
    "getEleLoadClassTags",
    "getEleLoadTags",
    "getEleLoadData",
    "getNodeLoadTags",
    "getNodeLoadData",
    "randomVariable",
    "getRVTags",
    "getRVParamTag",
    "getRVValue",
    "getMean",
    "getStdv",
    "getPDF",
    "getCDF",
    "getInverseCDF",
    "correlate",
    "performanceFunction",
    "gradPerformanceFunction",
    "transformUtoX",
    "wipeReliability",
    "updateMaterialStage",
    "sdfResponse",
    "probabilityTransformation",
    "startPoint",
    "randomNumberGenerator",
    "reliabilityConvergenceCheck",
    "searchDirection",
    "meritFunctionCheck",
    "stepSizeRule",
    "rootFinding",
    "functionEvaluator",
    "gradientEvaluator",
    "getNumThreads",
    "setNumThreads",
    "logFile",
    "setStartNodeTag",
    "hystereticBackbone",
    "stiffnessDegradation",
    "strengthDegradation",
    "strengthControl",
    "unloadingRule",
    "partition",
    "pressureConstraint",
    "domainCommitTag",
    "runFOSMAnalysis",
    "findDesignPoint",
    "runFORMAnalysis",
    "getLSFTags",
    "runImportanceSamplingAnalysis",
    "IGA",
    "NDTest",

    "tcl"
]

def _as_str_arg(arg):
    """
    Convert arg to a string that represents
    Tcl semantics.
    """
    if isinstance(arg, list):
        return f"{{{''.join(_as_str_arg(a) for a in arg)}}}"

    elif isinstance(arg, bool):
        return str(int(arg))

    elif isinstance(arg, dict):
        return "{\n" + "\n".join([
          f"{cmd} " + " ".join(_as_str_arg(a) for a in val)
              for cmd, val in kwds
    ]) + "}"
    else:
        return str(arg)



class OpenSeesPy:
    """
    This class is meant to be instantiated as a global singleton
    that is private to this Python module.

    It encapsulates an instance of Interpreter which implements an
    OpenSees state.
    """
    def __init__(self, *args, save=False, echo_file=None, **kwds):
        self._interp = Interpreter(*args,  **kwds)
        self._partial = partial
        self._save = save
        self._echo = echo_file # sys.stdout

        # Enable OpenSeesPy command behaviors
        self._interp.eval("pragma openseespy")

    def eval(self, cmd: str) -> str:
        "Evaluate a Tcl command"
        if self._echo is not None:
            print(cmd, file=self._echo)
        return self._interp.eval(cmd)

    def _str_call(self, proc_name: str, *args, **kwds)->str:
        """
        Invoke the Interpreter's eval method, calling
        a procedure named `proc_name` with arguments
        from args and kwds, after converting Python semantics
        to Tcl semantics (via _as_str_arg).

        For example, key-word arguments contained in the `kwds`
        dict are converted to a sequence of "-key" and "value"
        strings.
        """

        tcl_args = [_as_str_arg(i) for i in args]
        tcl_args += [
          f"-{key} " + _as_str_arg(val)
              for key, val in kwds.items()
        ]
        cmd = f"{proc_name} " + " ".join(tcl_args)

        # TODO: make sure errors print nicely
        try:
            ret = self.eval(cmd)
        except Exception as e:
            raise e

        if ret is None or ret == "":
            return None

        parts = ret.split()
        # Use json parse to cast return values from string. 
        # This is faster than the standard ast module.
        if len(parts) > 1:
            try:    return json.loads("[" + ",".join(parts) + "]")
            except: return ret
        else:
            try:    return json.loads(ret)
            except: return ret

    def pattern(self, *args, **kwds):
        self._current_pattern = args[1]
        return self._str_call("pattern", *args, **kwds)

    def load(self, *args, **kwds):
        pattern = self._current_pattern
        return self._str_call("nodalLoad", *args, "-pattern", pattern, **kwds)

    def section(self, *args, **kwds):
        self._current_section = args[1]
        # TODO: error handling
        return self._str_call("section", *args, **kwds)

    def patch(self, *args, **kwds):
        section = self._current_section
        return self._str_call("patch", *args, "-section", section, **kwds)

    def layer(self, *args, **kwds):
        section = self._current_section
        return self._str_call("layer", *args, "-section", section, **kwds)

    def fiber(self, *args, **kwds):
        section = self._current_section
        return self._str_call("fiber", *args, "-section", section, **kwds)


class Model:
    def __init__(self, *args, echo_file=None, **kwds):
        self._openseespy = OpenSeesPy(echo_file=echo_file)
        self._openseespy._str_call("model", *args, **kwds)

    def __getattr__(self, name: str):
        if name in {"pattern", "load", "eval", "section", "patch", "layer", "fiber"}:
            return getattr(self._openseespy, name)
        else:
            return self._openseespy._partial(self._openseespy._str_call, name)


# The global singleton, for backwards compatibility
_openseespy = OpenSeesPy()
# 
_tcl_echo   = None

def __getattr__(name: str):
    # For reference:
    #   https://peps.python.org/pep-0562/#id4
    if name in {"pattern", "load", "eval", "section", "patch", "layer", "fiber"}:
        return getattr(_openseespy, name)
    else:
        return _openseespy._partial(_openseespy._str_call, name)

