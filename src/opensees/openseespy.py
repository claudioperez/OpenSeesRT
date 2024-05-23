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

# something to compare the output of model.analyze to:
successful = 0

# A list of symbol names that are importable
# from this module. All of these are dynamically
# resolved by the function __getattr__ below.
__all__ = [
# 
    "tcl",
    "OpenSeesError",

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

_PROTOTYPES = {
}

# Commands that are pre-processed in Python
# before forwarding to the Tcl interpreter
_OVERWRITTEN = {
    "timeSeries",
    "pattern", "load",
    "eval",
    "section", "patch", "layer", "fiber",
    "block2D",
    "mesh"
}


class OpenSeesError(Exception):
    pass

def _split_iter(source, sep=None, regex=False):
    """
    generator version of str.split()

    :param source:
        source string (unicode or bytes)

    :param sep:
        separator to split on.

    :param regex:
        if True, will treat sep as regular expression.

    :returns:
        generator yielding elements of string.
    """
    if sep is None:
        # mimic default python behavior
        source = source.strip()
        sep = "\\s+"
        if isinstance(source, bytes):
            sep = sep.encode("ascii")
        regex = True

    if regex:
        # version using re.finditer()
        if not hasattr(sep, "finditer"):
            sep = re.compile(sep)
        start = 0
        for m in sep.finditer(source):
            idx = m.start()
            assert idx >= start
            yield source[start:idx]
            start = m.end()
        yield source[start:]

    else:
        # version using str.find(), less overhead than re.finditer()
        sepsize = len(sep)
        start = 0
        while True:
            idx = source.find(sep, start)
            if idx == -1:
                yield source[start:]
                return
            yield source[start:idx]
            start = idx + sepsize


class OpenSeesPy:
    """
    This class is meant to be instantiated as a global singleton
    that is private to this Python module.

    It encapsulates an instance of Interpreter which implements an
    OpenSees state.
    """
    def __init__(self, *args, save=False, echo_file=None, **kwds):
        self._interp  = Interpreter(*args,  **kwds)
        self._partial = partial
        self._save    = save
        self._echo    = echo_file

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
            try:    return list(map(json.loads, parts)) #json.loads("[" + ",".join(parts) + "]")
#           try:    return json.loads("[" + ",".join(parts) + "]")
            except: return ret

        elif proc_name == "eigen":
            # "eigen" should always return a list
            return [float(ret)]

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
        """
        ['Path', 1, '-values', 0.0, 5.0, 8.0, 7.0, 5.0, 3.0, 2.0, 1.0, 0.0, '-time', 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        ['Path', 1, '-values', [0.0, 5.0, 8.0, 7.0, 5.0, 3.0, 2.0, 1.0, 0.0], '-time', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]]
        """

        args = list(args)
        if "-values" in args:
            iv = args.index("-values")
            # Count the number of floating-point arguments
            for nv, value in enumerate(args[iv+1:]):
                if not isinstance(value, float):
                    nv += 1
                    break
            else:
                # if we didnt break out of the for loop
                nv += 2

            values = args[iv+1:iv+nv]
            args = [a for a in args[:iv+1]] + [values] + [a for a in args[iv+nv:]]

        if "-time" in args:
            it = args.index("-time")
            for nt, value in enumerate(args[it+1:]):
                if not isinstance(value, float):
                    nt += 1
                    break
            else:
                # if we didnt break out of the for loop
                nt += 2

            time = args[it+1:it+nt]
            args = [a for a in args[:it+1]] + [time] + [a for a in args[it+nt:]]

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

    def mesh(self, type, tag: int, *args, **kwds):
        if type == "line":
            return self._mesh_line(tag, 2, args[1:3], *args[3:7], args[7:])

    def _mesh_line(self, tag, numnodes, ndtags, id, ndf:int, meshsize, eleType='', eleArgs=()):
        import numpy as np
        from itertools import count

        ndI, ndJ = ndtags
        add_node    = partial(self._str_call, "node")
        add_element = partial(self._str_call, "element")

        xi = np.array(self._str_call("nodeCoord", ndI))
        xj = np.array(self._str_call("nodeCoord", ndJ))

        L  = np.linalg.norm(xj - xi)
        nn = int(L//meshsize) + 1

        nodes = [None for _ in range(nn)]
        nodes[0]    = ndI
        nodes[nn-1] = ndJ

        node_tags = set(self._str_call("getNodeTags"))
        new_node  = filter(lambda i: i not in node_tags, count(1))
        elem_tags = set(self._str_call("getEleTags") or [])
        new_elem  = filter(lambda i: i not in elem_tags, count(1))

        for i,x in enumerate(np.linspace(xi, xj, nn, endpoint=True)[1:]):

            node_tag = next(new_node)
            add_node(node_tag, *x)

            nodes[i+1] = node_tag

            elem_tag = next(new_elem)

            if i < nn:
                add_element(eleType,elem_tag,nodes[i],nodes[i+1],*eleArgs)



    def section(self, type: str, sec_tag: int, *args, **kwds):
        self._current_section = sec_tag
        # TODO: error handling

        if "shape" in kwds:
            from opensees.section import from_shape
            ndm = int(self.eval("getNDM"))
            # kwds["shape"] looks like ("W14X90", matTag, (20,4), units?)
            shape = from_shape(type, *kwds.pop("shape"), ndm=ndm)
        else:
            shape = None

        ret = self._str_call("section", type, sec_tag, *args, **kwds)

        if shape is not None:
            for fiber in shape.fibers:
                self._str_call("fiber", *fiber.coord, fiber.area, fiber.material, section=sec_tag)

        return ret

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

    def export(self, *args, **kwds):
        return self._openseespy._interp.export(*args, **kwds)

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
        return "{\n" + "\n".join([
          f"{cmd} " + " ".join(_as_str_arg(a) for a in val)
              for cmd, val in arg.items()
        ]) + "}"

    else:
        return str(arg)


class FedeasModel(Model):
    @property
    def nf(self):
        return self.numDOF()

# The global singleton, for backwards compatibility
_openseespy = OpenSeesPy()


def __getattr__(name: str):
    # For reference:
    #   https://peps.python.org/pep-0562/#id4
    if name in _OVERWRITTEN:
        return getattr(_openseespy, name)
    else:
        return _openseespy._partial(_openseespy._str_call, name)

