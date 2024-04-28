"""
This module implements the OpenSeesPy interface.
Imports can be performed exactly as one would
from openseespy, for example:

>>> import opensees.openseespy as ops

>>> from opensees.openseespy import node, model

>>> from opensees.openseespy import *

"""
import json
from functools import partial

from .tcl import Interpreter

# A list of symbol names that are importable
# from this module. All of these are dynamically
# resolved by the function __getattr__ below.
__all__ = [
# 
    "tcl"
    "OpenSeesError"

# OpenSeesPy attributes
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
#   "runFOSMAnalysis",
    "findDesignPoint",
    "runFORMAnalysis",
    "getLSFTags",
    "runImportanceSamplingAnalysis",
    "IGA",
    "NDTest",
]

# 
_OVERWRITTEN = {
    "timeSeries",
    "pattern", "load",
    "eval",
    "section", "patch", "layer", "fiber",
    "block2D"
}


class OpenSeesError(Exception):
    pass


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

    def _str_call(self, proc_name: str, *args, _final=None, **kwds)->str:
        """
        Invoke the Interpreter's eval method, calling
        a procedure named `proc_name` with arguments
        from args and kwds, after converting Python semantics
        to Tcl semantics (via _as_str_arg).

        For example, key-word arguments contained in the `kwds`
        dict are converted to a sequence of "-key" and "value"
        strings.
        """

        tcl_args = (_as_str_arg(i) for i in args)
        tcl_kwds = (
          (f"-{key.replace('_','-')}" if val else "") if isinstance(val, bool)
          else f"-{key} " + _as_str_arg(val)
              for key, val in kwds.items()
        )
        cmd = f"{proc_name} {' '.join(tcl_args)} {' '.join(tcl_kwds)}"

        if _final is not None:
            cmd += _as_str_arg(_final)

        # TODO: make sure errors print nicely
        try:
            ret = self.eval(cmd)
        except Exception as e:
            raise OpenSeesError() from e

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


    def eval(self, cmd: str) -> str:
        "Evaluate a Tcl command"
        if self._echo is not None:
            print(cmd, file=self._echo)
        return self._interp.eval(cmd)


    def block2D(self, *args, **kwds):
        if isinstance(args[5], list):
            return self._str_call("block2D", *args, **kwds)

        # We have to imitate the OpenSeesPy parser, which
        # *requires* hard-coding the number of element args
        # expected by each element type. This is terribly
        # unstable and limited and should only be used when 
        # backwards compatibility with the original OpenSeesPy 
        # is absolutely necessary.
        elem_name = args[4]
        elem_argc = {
            "quad":         9,
            "stdquad":      9,

            "shell":        7,
            "shellmitc4":   7,

            "shellnldkgq":  7,
            "shelldkgq":    7,

            "bbarquad":     8,

            "enhancedquad": 9,

            "sspquad":      9
        }[elem_name.lower()] -1

        elem_args = list(args[5:elem_argc])

        nl  = '\n'
        ndm = self._str_call("getNDM")
        # loop over remaining args to form node coords
        node_args = f"""{{
            {nl.join(" ".join(map(str,args[elem_argc+i*(ndm+1):elem_argc+(i+1)*(ndm+1)])) for i in range(int(len(args[elem_argc:])/(ndm+1))))}
        }}"""

        return self._str_call("block2D", *args[:5], elem_args, node_args)


    def timeSeries(self, *args, **kwds):
        if "-values" in args:
            iv = list(args).index("-values")
            values = list(args[iv+1:])
            args = [a for a in args[:iv+1]] + [values]

        return self._str_call("timeSeries", *args, **kwds)

    def pattern(self, *args, load=None, **kwds):
        self._current_pattern = args[1]

        if load is not None:
            loads = [
                    ("load", k, *v) for k,v in load.items()
            ]
            return self._str_call("pattern", *args, **kwds, _final=loads)
        else:
            return self._str_call("pattern", *args, **kwds)

    def load(self, *args, pattern=None, load=None, **kwds):
        if pattern is None:
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

    def asdict(self):
        """April 2024"""
        return self._openseespy._interp.serialize()

    def setFactor(self, factor):
        pass

    def getIterationCount(self):
        return self._openseespy._str_call("numIter")

    def getResidual(self):
        return self._openseespy._str_call("printB", "-ret")

    def getTangent(self, **kwds):
        import numpy as np
        A = np.array(self._openseespy._str_call("printA", "-ret", **kwds))
        return A.reshape([int(np.sqrt(len(A)))]*2)

    def __getattr__(self, name: str):
        if name in _OVERWRITTEN:
            return getattr(self._openseespy, name)
        else:
            return self._openseespy._partial(self._openseespy._str_call, name)

def _as_str_arg(arg, name: str = None):
    """
    Convert arg to a string that represents
    Tcl semantics.
    """
    import numpy as np
    if isinstance(arg, (list,np.ndarray)):
        return f"{{{' '.join(_as_str_arg(a) for a in arg)}}}"

    elif isinstance(arg, tuple):
        return " ".join(map(str, arg))

    # parse commands like `section Fiber {...}`
    elif isinstance(arg, dict):
        print(arg)
        return "{\n" + "\n".join([
          f"{cmd} " + " ".join(_as_str_arg(a) for a in val)
              for cmd, val in arg.items()
        ]) + "}"

    else:
        return str(arg)




# The global singleton, for backwards compatibility
_openseespy = OpenSeesPy()


def __getattr__(name: str):
    # For reference:
    #   https://peps.python.org/pep-0562/#id4
    if name in _OVERWRITTEN:
        return getattr(_openseespy, name)
    else:
        return _openseespy._partial(_openseespy._str_call, name)

