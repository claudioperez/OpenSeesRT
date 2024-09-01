"""
Matrices. (:mod:`matrices`)
=====================================================

.. currentmodule:: matrices

This module provides functions and classes for constructing various
structural analysis matrices. Matrices are generally constructed with 
functions that are defined outside of classes (as opposed to using the 
__init__ method of a class) in an effort to present basic structural analysis 
concepts with minimal clutter from computer science concepts which may
be unfamiliar to engineers.

"""

import numpy as np
import pandas as pd
import scipy.linalg
import matplotlib.pyplot as plt


import numpy as np
import matplotlib.pyplot as plt
import copy


def _create_model(model):
    pass

def static_matrix(model):
    pass

def strain_matrix(model):
    pass

##<****************************Model objects************************************
# Model objects/classes
# These should be created using the methods above
#******************************************************************************
class _Node:
    def __init__(self, model, tag: str, ndf, xyz):

        self.xyz = np.array([xi for xi in xyz if xi is not None])

        self.tag = tag
        self.xyz0 = self.xyz # coordinates in base configuration (unstrained, not necessarily unstressed).
        self.xyzi = self.xyz # coordinates in reference configuration.  

        self.x0: float = xyz[0] # x-coordinate in base configuration (unstrained, not necessarily unstressed).  
        self.y0: float = xyz[1] # y-coordinate in base configuration (unstrained, not necessarily unstressed).  
        self.z0: float = xyz[2] # z-coordinate in base configuration (unstrained, not necessarily unstressed).  


        self.x: float = xyz[0]
        self.y: float = xyz[1]
        self.z: float = xyz[2]


        self.rxns = [0]*ndf
        self.model = model
        self.elems = []

        self.p = {dof:0.0 for dof in model.ddof}

    def __repr__(self):
        return 'nd-{}'.format(self.tag)

    def p_vector(self):
        return np.array(list(self.p.values()))

    @property
    def dofs(self):
        # if self.model.DOF == None: self.model.numDOF()
        idx = self.model.nodes.index(self)
        return np.array(self.model.DOF[idx],dtype=int)


class Rxn:
    def __init__(self, node, dirn):
        self.node = node
        self.dirn = dirn

    def __repr__(self):
        return 'rxn-{}'.format(self.dirn)


class Hinge:
    def __init__(self, elem, node):
        self.elem = elem
        self.node = node


class Model:
    def __init__(self, ndm, ndf):
        """Basic model object

        Parameters
        -----------
        ndm:
            number of model dimensions
        ndf:
            number of degrees of freedom (dofs) at each node

        """
        self.ndf: int = ndf
        self.ndm: int = ndm
        self.DOF: list = None

        # Define DOF list indexing 
        if ndm == 1:
            self.prob_type = '1d'
            self.ddof: dict = {'x': 0}  # Degrees of freedom at each node

        if ndf == 2:
            self.prob_type = '2d-truss'
            self.ddof: dict = { 'x': 0, 'y': 1} # Degrees of freedom
        elif ndm == 2 and ndf ==3:
            self.prob_type = '2d-frame'
            self.ddof: dict = { 'x': 0, 'y': 1, 'rz':2}
        elif ndm == 3 and ndf ==3:
            self.prob_type = '3d-truss'
            self.ddof: dict = { 'x': 0, 'y': 1, 'z':2}
        elif ndm == 3 and ndf ==6:
            self.prob_type = '3d-frame'
            self.ddof: dict = { 'x': 0, 'y': 1, 'z':2, 'rx':3, 'ry':4, 'rz':5}

        # model inventory lists
        self.elems: list = []
        self.nodes: list = []
        self.rxns:  list = []
        self.hinges: list = []
        self.iforces: list = []
        self.loads: list = []
        self.states: list = []
        self.redundants: list = []

        # model inventory dictionaries
        self.delems: dict = {}
        self.dnodes: dict = {}
        self.dxsecs: dict = {}
        self.dhinges: dict = {}
        self.dstates: dict = {}
        self.dxsecs: dict = {}
        self.materials: dict = {}
        self.xsecs: dict = {}
        self.dredundants: dict = {}

        # Initialize default material/section properties
        self.material('default', 1.0)
        self.xsection('default', 1.0, 1.0)

    @property
    def rel(self):
        return [rel for elem in self.elems for rel in elem.rel.values()]

    @property
    def nn(self) -> int:
        """return number of nodes in model"""
        return len(self.nodes)

    @property
    def nr(self) -> int:
        """return number of constrained dofs in model"""
        return len(self.rxns)

    @property
    def ne(self) -> int:
        """return number of elements in model"""
        return len(self.elems)

    @property
    def nQ(self):
        f = 0
        for elem in self.elems:
            f += len(elem.basic_forces)
        return f

    @property
    def nq(self):
        f = []
        for elem in self.elems:
            f.append(sum([1 for x in elem.q]))
        return f

    @property
    def nv(self):
        lst = []
        for elem in self.elems:
            lst.append(sum([1 for x in elem.v]))
        return lst

    @property
    def nf(self) -> int:
        x = self.nt - self.nr
        return x

    @property
    def nt(self) -> int:
        return self.ndf*self.nn

    @property
    def rdofs(self): 
        """Return list of restrained dofs"""
        DOF = self.DOF

        return []

    @property
    def NOS(self) -> int:
        nf = self.nf
        nq = sum(self.nq)
        return nq - nf

    @property
    def basic_forces(self):
        return np.array([q for elem in self.elems for q in elem.basic_forces ])

    @property
    def rdnt_forces(self):
        cforces = self.cforces
        return np.array([q  for q in cforces if q.redundant ])

    @property
    def cforces(self):
        return np.array([q for elem in self.elems for q in elem.basic_forces if not q.rel])

    @property 
    def eforces(self):
        """Return array of elastic element forces"""
        return np.array([q for elem in self.elems for q in elem.basic_forces if (q.plastic_event is None)])

    @property
    def idx_c(self):
        cforces = self.cforces
        forces = self.basic_forces
        idx_c = np.where(np.isin(forces,cforces))[0]
        return idx_c

    @property
    def idx_e(self):
        """return indices of elastic basic (not plastic) forces"""
        cforces = self.cforces
        eforces = self.eforces
        idx_e = np.where(np.isin(cforces,eforces))[0]
        return idx_e

    @property
    def idx_f(self):
        return np.arange(0,self.nf)
    @property
    def idxx_f(self):
        return None

    @property
    def idx_i(self):
        rdts = self.rdnt_forces
        forces = self.basic_forces
        idx_i = np.where(np.logical_not(np.isin(forces, rdts)))[0]
        return idx_i

    @property
    def idx_x(self):
        rdts = self.rdnt_forces
        forces = self.basic_forces
        idx_x = np.where(np.isin(forces, rdts))[0]
        return idx_x

    def node(self, tag: str, x: float, y=None, z=None, mass: float=None):
        newNode = _Node(self, tag, self.ndf, [x, y, z], mass)
        self.nodes.append(newNode)
        self.dnodes.update({newNode.tag : newNode})
        return newNode

    def state(self, method="Linear"):
        if self.DOF==None: self.numDOF()
        newState = State(self, method)

        ElemTypes = {type(elem) for elem in self.elems}
        StateVars = {key for elem in ElemTypes for key in elem.stateVars.keys() }

        stateDict = {
                var : {
                        elem.tag : copy.deepcopy(elem.stateVars[var])
                            for elem in self.elems if var in elem.stateVars.keys()
                } for var in StateVars
        }
        self.states.append(stateDict)
        return stateDict

    def numDOF(self):
        crxns = self.ndf*len(self.nodes) - len(self.rxns)+1
        df = 1
        temp = []
        for node in self.nodes:
            DOFs = []
            for rxn in node.rxns:
                if not(rxn):
                    DOFs.append(df)
                    df += 1
                else:
                    DOFs.append(crxns)
                    crxns += 1
            temp.append(DOFs)
        self.DOF = temp
        return self.DOF



    def fix(self, node, dirn=["x","y","rz"]): # for dirn enter string (e.g. "x", 'y', 'rz')
        if isinstance(dirn,list):
            rxns = []
            for df in dirn:
                newRxn = Rxn(node, df)
                self.rxns.append(newRxn)
                rxns.append(newRxn)
                node.rxns[self.ddof[df]] = 1
            return rxns
        else:
            newRxn = Rxn(node, dirn)
            self.rxns.append(newRxn)
            node.rxns[self.ddof[dirn]] = 1
            return newRxn

 # Other
    def material(self, tag: str, E: float):
        newMat = Material(tag, E)
        self.materials[tag]=newMat
        return newMat

    def xsection(self, tag: str, A: float, I: float):
        newXSect = XSect(tag, A, I)
        self.xsecs[tag] = newXSect
        return newXSect

 # Elements
    def add_element(self, element):
        """Add a general element to model

        Parameters
        ---------
        element : obj

        """

        self.elems.append(element)
        self.delems.update({element.tag:element})

        for node in element.nodes:
            node.elems.append(element)

        return element

    def add_elements(self, elements):
        """Add a general element to model

        Parameters
        ---------
        element : obj

        """
        for element in elements:
            self.elems.append(element)
            self.delems.update({element.tag:element})
            for node in element.nodes:
                node.elems.append(element)

        return elements


    def beam(self, tag: str, iNode, jNode, mat=None, sec=None, Qpl=None,**kwds):
        """Add a 2D linear Euler-Bernouli beam object to model

        Parameters
        ---------
        tag : str
            string used for identifying object
        iNode : ema.Node
            node object at element i-end
        jNode : ema.Node
            node object at element j-end
        mat : ema.Material 

        sec : ema.Section


        """

        if mat is None:
            E = kwds["E"] if "E" in kwds else self.materials["default"].E
        else:
            E = mat.E

        if sec is None:
            A = kwds["A"] if "A" in kwds else self.xsecs["default"].A
            I = kwds["I"] if "I" in kwds else self.xsecs["default"].I
        else:
            A = sec.A
            I = sec.I

        from ema.elements import Beam

        newElem = Beam(tag, iNode, jNode, E, A, I)
        self.elems.append(newElem)
        self.delems.update({newElem.tag:newElem})
        # self.connect([iNode, jNode], "Beam") # considering deprecation
        iNode.elems.append(newElem)
        jNode.elems.append(newElem)

        if Qpl is not None:
            newElem.Qpl = np.zeros((3,2))
            newElem.Np = [Qpl[0], Qpl[0]]
            newElem.Mp = [[Qpl[1], Qpl[1]],[Qpl[2], Qpl[2]]]
            for i, key in enumerate(newElem.Qp['+']):
                newElem.Qp['+'][key] = newElem.Qp['-'][key] = Qpl[i] # consider depraction of elem.Qp in favor of elem.Qpl
                newElem.Qp['+'][key] = newElem.Qp['-'][key] = Qpl[i]
                newElem.Qpl[i,:] = Qpl[i] # <- consider shifting to this format for storing plastic capacities
        return newElem


    def truss(self, tag: str, iNode, jNode, mat=None, xsec=None, Qpl=None,A=None,E=None):
        from ema.elements import Truss

        if mat is None: mat = self.materials['default']
        if E is None: E = mat.E
        # cross section
        if xsec is None: xsec = self.xsecs['default']
        if A is None: A = xsec.A


        newElem = Truss(tag, iNode, jNode, E, A)
        self.delems.update({newElem.tag:newElem})
        self.elems.append(newElem)
        iNode.elems.append(newElem)
        jNode.elems.append(newElem)

        if Qpl is not None:
            newElem.Np = [Qpl[0], Qpl[0]]
            newElem.Qp['+']['1'] = newElem.Qp['-']['1'] = Qpl[0]
        return newElem


    def taprod(self, tag: str, iNode, jNode, mat=None, xsec=None, Qpl=None,A=None,E=None):
        """Construct a tapered rod element with variable E and A values."""
        from ema.elements import TaperedTruss

        if mat is None:
            mat = self.materials['default']
        if E is None:
            E = mat.E
        # cross section
        if xsec is None: xsec = self.xsecs['default']
        if A is None: A = xsec.A

        newElem = TaperedTruss(tag, iNode, jNode, E, A)
        self.delems.update({newElem.tag:newElem})
        self.elems.append(newElem)
        iNode.elems.append(newElem)
        jNode.elems.append(newElem)

        if Qpl is not None:
            newElem.Np = [Qpl[0], Qpl[0]]
            newElem.Qp['+']['1'] = newElem.Qp['-']['1'] = Qpl[0]
        return newElem

    def truss3d(self, tag: str, iNode, jNode, mat=None, xsec=None):
        """Add an ema.Truss3d object to model

        Parameters
        ---------

        """
        from ema.elements import Truss3D
        if mat is None: mat = self.materials['default']
        if xsec is None: xsec = self.xsecs['default']
        newElem = Truss3D(tag, iNode, jNode, mat, xsec)
        self.elems.append(newElem)
        self.delems.update({newElem.tag:newElem})
        iNode.elems.append(newElem)
        jNode.elems.append(newElem)
        return newElem


    def hinge(self, elem, node):
        "pin a beam end."
        newHinge = Hinge(elem, node)
        self.hinges.append(newHinge)
        if node == elem.nodes[0]:
            elem.rel['2'] = True
            elem.q.pop('2')
            elem.basic_forces[1].rel = True
        elif node == elem.nodes[1]:
            elem.rel['3'] = True
            elem.q.pop('3')
            elem.basic_forces[2].rel = True
        else:
            raise Exception("element {} is not bound to node {}.".format(elem, node))

        return newHinge


    def redundant(self, elem, nature):
        newq = IntForce(elem, nature)
        elem.red[nature] = True
        self.redundants.append(newq)


class rModel(Model):
    def __init__(self, ndm, ndf):
        super().__init__(ndm=2, ndf=3)
        self.material('default', 1.0)
        self.xsection('default', 1.0, 1.0)

    def isortho(self, elem):
        if (abs(elem.cs) == 1.0) or (abs(elem.sn) == 1.0):
            return True
        else:
            return False

    def numdofs(self):
        current_rxn = 1
        current_dof = 1
        rxn_ixs = []
        DOFs = [[0, 0, 0] for node in self.nodes]
        for i, node in enumerate(self.nodes):
            # x-dof
            dirn = 0
            if not(node.rxns[dirn]): # node is free
                if not(DOFs[i][dirn]): # node unassigned
                    for elem in node.elems:
                        if abs(elem.cs) == 1.0: # x-dof coupled to far end
                            if elem.nodes[0] == node:
                                far_node = self.nodes.index(elem.nodes[1])
                            if elem.nodes[1] == node:
                                far_node = self.nodes.index(elem.nodes[0])

                            if not(DOFs[far_node][dirn]): # Far node dof unassigned
                                if not(self.nodes[far_node].rxns[dirn]): # Far node is free
                                    DOFs[far_node][dirn] = current_dof
                                    DOFs[i][dirn] = current_dof
                                    current_dof += 1
                                else: # Far node is fixed
                                    DOFs[far_node][dirn] = current_rxn
                                    DOFs[i][dirn] = current_rxn
                                    current_rxn += 1
                                    rxn_ixs.append( (i,dirn) )
                                    rxn_ixs.append( (far_node,dirn) )
                            else: # Far node dof already assigned
                                DOFs[i][dirn] = DOFs[far_node][dirn]
                                if self.nodes[far_node].rxns[dirn]: # Far node is fixed
                                    rxn_ixs.append( (i,dirn) )

                        elif all([abs(elem.cs) != 1.0 for elem in node.elems]): # x-dof free/uncoupled
                            if not DOFs[i][dirn]:
                                DOFs[i][dirn] = current_dof
                                current_dof += 1

            else: # node is fixed
                if not DOFs[i][dirn]: # node is unassigned
                    DOFs[i][dirn] = current_rxn 
                    current_rxn += 1
                    rxn_ixs.append( (i,dirn) )

            # y-dof
            dirn = 1
            if not(node.rxns[dirn]):
                if not(DOFs[i][dirn]):
                    for elem in node.elems:
                        if abs(elem.sn) == 1.0:
                            # get far node index
                            if elem.nodes[0] == node:
                                far_node = self.nodes.index(elem.nodes[1])
                            if elem.nodes[1] == node:
                                far_node = self.nodes.index(elem.nodes[0])

                            if not(DOFs[far_node][dirn]):
                                if not(self.nodes[far_node].rxns[dirn]):
                                    DOFs[far_node][dirn] = current_dof
                                    DOFs[i][dirn] = current_dof
                                    current_dof += 1
                                else:
                                    DOFs[far_node][dirn] = current_rxn
                                    DOFs[i][dirn] = current_rxn
                                    current_rxn += 1
                                    rxn_ixs.append( (i,dirn) )
                                    rxn_ixs.append( (far_node,dirn) )
                            else: 
                                DOFs[i][dirn] = DOFs[far_node][dirn]
                                if self.nodes[far_node].rxns[dirn]:
                                    rxn_ixs.append( (i,dirn) )
                        elif all([abs(elem.sn) != 1.0 for elem in node.elems]):
                            if not(DOFs[i][dirn]):
                                DOFs[i][dirn] = current_dof
                                current_dof += 1
            else:
                if not(DOFs[i][dirn]):
                    DOFs[i][dirn] = current_rxn
                    current_rxn += 1
                    rxn_ixs.append( (i,dirn) )

          # rz-dof
            dirn = 2
            if not(node.rxns[2]):
                DOFs[i][dirn] = current_dof
                current_dof += 1
            else:
                DOFs[i][dirn] = current_rxn
                current_rxn += 1
                rxn_ixs.append( (i,dirn) )

        for ids in rxn_ixs:
            DOFs[ids[0]][ids[1]] += current_dof - 1

        return DOFs


    def numDOF(self):
        crxns = self.ndf*len(self.nodes) - len(self.rxns)+1
        df = 1
        temp = []
        for node in self.nodes:
            DOFs = []
            for rxn in node.rxns:
                if not(rxn):
                    DOFs.append(df)
                    df += 1
                else:
                    DOFs.append(crxns)
                    crxns += 1
            temp.append(DOFs)
        self.DOF = temp
        return self.DOF

    @property
    def triv_forces(self):
        """list of trivial axial forces"""
        lst = []
        for elem in self.elems:
            if len(elem.basic_forces) > 1:
                if elem.dofs[0]==elem.dofs[3] or elem.dofs[1]==elem.dofs[4]:
                    lst.append(elem.basic_forces[0])
        return np.array(lst)

    # @property
    # def basic_forces(self):
    #     # bmax = self.TrAx_forces
    #     forces = np.array([q for elem in self.elems for q in elem.basic_forces if not q.rel])
    #     return forces

    @property
    def cforces(self):
        triv = self.triv_forces
        arry = np.array([q for elem in self.elems for q in elem.basic_forces if (
            (q.plastic_event is None) and (
            not q in triv) and (
            not q.rel))])
        return arry

    @property
    def nr(self):
        return len(self.rxns)

    # @property
    # def nq(self):
    #     f = []
    #     for elem in self.elems:
    #         f.append(sum([1 for q in elem.basic_forces if not q.rel and (not q in self.triv_forces)]))
    #     return f

    @property
    def nv(self):
        """Returns number of element deformations in model"""
        lst = []
        for elem in self.elems:
            lst.append(sum([1 for x in elem.v]))
        return lst



    @property
    def fdof(self): 
        """Return list of free dofs"""
        pass

    @property
    def nt(self):
        nt = max([max(dof) for dof in self.DOF])
        return nt

    @property
    def nm(self):
        """No. of kinematic mechanisms, or the no. of dimensions spanned by the 
        m independent inextensional mechanisms.

        """
        # A = A_matrix(mdl)
        pass

    @property
    def NOS(self):
        nf = self.nf
        nq = sum(self.nq)
        return nq - nf



class Material():
    def __init__(self, tag, E, nu=None):
        self.tag = tag
        self.E: float = E
        self.elastic_modulus: float = E 
        self.poisson_ratio = nu

class XSect():
    def __init__(self, tag, A, I):
        self.tag = tag
        self.A: float = A
        self.I: float = I

def  UnitProperties(): # Legacy; consider removing and adding unit materials to model by default
    return (Material('unit', 1.0), XSect('unit', 1.0, 1.0))

class State:
    """STATE is a data structure with information about the current state of the structure in fields

    """
    def __init__(self, model, method="Linear"):
        self.model = model
        self.num = len(model.states)
        self.data = {"Q": [[qi for qi in elem.q.values()] for elem in model.elems],
                     "P": {str(dof):0 for dof in [item for sublist in model.DOF for item in sublist]},
                     "DOF": 'model.numDOF(model)'
                    }
        self.method = method


    def eload(self, elem, mag, dirn='y'):
        if type(elem) is str:
            elem = self.model.delems[elem]
        if not(type(mag) is list):
            if dirn=='y':
                mag = [0.0, mag]
            elif dirn=='x':
                mag = [mag, 0.0]
        elem.w[self.num] = mag

##>*****************************************************************************


#
# ##>*****************************************************************************
#


settings = {
    "DATAFRAME_LATEX": True, 
}

subscripts = {
    'part': 'subscript',
     None: '',
    'initial': '0',
    'continuous':'c',
    'primary':'i',
    'reactions':'d',
    'elem-load':'w',
    'free': 'f',
}

def del_zeros(mat):
    delrows = np.where(~mat.any(axis=1))[0]
    delcols = np.where(~mat.any(axis=0))[0]

    newM = np.delete(mat, delcols, axis=1)
    newM = np.delete(newM, delrows, axis=0)
    return newM

def _elem_dofs(Elem):
    dofs = []
    for node in Elem.nodes:
        dofs.extend(node.dofs)
    return dofs

def transfer_vars(item1, item2):
    for key, value in item1.__dict__.items():
        item2.__dict__[key] = value

def Localize(U_vector, P_vector, model=None):
    if model is None:
        model = U_vector.model
    A =  A_matrix(model)
    Q0 = Q0_vector(model)
    Ks = Ks_matrix(model)

    V = A.f @ U_vector.f
    Q  = Ks@V + Q0
    return V, Q



class Structural_Vector(np.ndarray):
    column_data = ["vector"]
    row_data = None
    subs = [None]
    tag = 'Vector'
    def __new__(cls, mat):
        mat = mat
        return np.asarray(mat).view(cls)

    def __init__(self, mat):
        self.tag = None

    def _repr_html_(self):
        try:
            out = self.df.to_html()
            return out
        except:
            return self

    def __add__(self, other):
        if isinstance(other, type(self)):
            out = np.add(self, other).view(type(self))
            transfer_vars(self, out)
        else:
            out = super().__add__(other)
        return out

    def get(self, key):
        idx = np.array([i for i,j in enumerate(self.row_data) if str(j) == key], dtype=int)
        # row = self.row_data.index(component)
        return self[idx]

    def set_item(self, key, value):
        idx = np.array([i for i,j in enumerate(self.row_data) if str(j) == key], dtype=int)
        self[idx] = value

    def rows(self, component):
        idx = np.where(np.isin(self.row_data,component))[0]
        newV = self[idx]
        newV.row_data = np.array(self.row_data)[idx]
        newV.model = self.model
        return newV

    @property
    def df(self):
        row_data = ['$'+str(tag)+'$' for tag in self.row_data]
        header = '$'+self.tag+'_{{'
        try:
            for sub in self.subs:
                header += subscripts[sub]
        except: pass
        print
        header += '}}$'

        return pd.DataFrame(np.around(self,14), index=row_data, columns=[header])

    @property
    def symb(self):
        var = []
        for eid in self.row_data:
            var.append(sp.symbols(self.tag+eid))
        return sp.Matrix(var)

    @property
    def disp(self):
        return sp.Matrix(self)


class Structural_Matrix(np.ndarray):
    column_data = None
    row_data = None
    c_ridx = None # row indexes for use in .c method
    c_cidx = None # column indexes for use in .c method
    tag = None

    def __new__(cls, mat):
        mat = mat
        return np.asarray(mat).view(cls)

    def __init__(self, mat):
        self.tag = None

    def __matmul__(self, other):
        if isinstance(other, Structural_Matrix):
            out = np.matmul(self,other).view(Structural_Matrix)
            transfer_vars(self,out)
            out.column_data = other.column_data

        elif isinstance(other, Structural_Vector):
            out = np.matmul(self,other).view(Structural_Vector)
            out.row_data = self.row_data
            # out.column_data = ['.']
        else:
            out = np.matmul(self,other).view(Structural_Vector)
            out.row_data = self.row_data
            # out.column_data = ['.']
        return out

    def __add__(self, other):
        out = np.add(self,other)
        if (isinstance(other, Structural_Matrix)
            or isinstance(other, Structural_Vector)):
            out.row_data = self.row_data
        return out



    def __truediv__(self, other):
        out = np.ndarray.__truediv__(self, other) 
        if (isinstance(other, float) 
            or isinstance(other, Structural_Matrix)
            or isinstance(other, Structural_Vector)):
            out = np.ndarray.__truediv__(self, other).view(Structural_Matrix)
            transfer_vars(self, out)
        else:
            out = np.ndarray.__truediv__(self, other)
        return out

    def _repr_html_(self):
        try: 
            df = self.df
            return df.to_html()
        except:
            try:
                return pd.DataFrame(self).to_html()
            except:
                pass

    @property
    def disp(self):
        return sp.Matrix(self)

    @property
    def df(self):
        if settings['DATAFRAME_LATEX']:
            row_data = ['$'+str(tag)+'$' for tag in self.row_data]
            column_data = ['$'+str(tag)+'$' for tag in self.column_data]
        else:
            row_data = [str(i) for i in self.row_data]
            column_data = [str(i) for i in self.column_data]
        return pd.DataFrame(np.around(self,5), index=row_data, columns=column_data)

    @property
    def inv(self):
        mat = np.linalg.inv(self)
        transfer_vars(self, mat)
        mat.row_data = self.column_data
        mat.column_data = self.row_data
        return mat

    @property
    def rank(self):
        """Return the rank of a matrix"""
        return np.linalg.matrix_rank(self)

    @property
    def lns(self):
        """Return a basis for the left nullspace of a matrix."""
        return scipy.linalg.null_space(self.T)

    @property
    def nls(self):
        """return a basis for the nullspace of matrix."""
        return scipy.linalg.null_space(self)

    @property
    def ker(self):
        "Return a basis for the kernel (nullspace) of a matrix."
        kernel = scipy.linalg.null_space(self) 
        ker = Structural_Matrix(kernel)
        transfer_vars(self,ker)
        ker.row_data = self.column_data 
        ker.column_data = [str(i+1) for i in range(len(ker[0]))]
        return ker

    def lu(self):

        return scipy.linalg.lu(self)

    @property
    def c(self):
        delcols = self.c_cidx
        delrows = self.c_ridx

        newM = np.delete(self, delrows, axis=0).view(type(self))
        if delcols: newM = np.delete(newM, delcols, axis=1).view(type(self))

        transfer_vars(self, newM)

        if delcols: newM.column_data = list(np.delete(self.column_data, delcols))

        newM.row_data = list(np.delete(self.row_data, delrows))
        return newM

    def round(self, num):
        newM = np.around(self, num).view(Structural_Matrix)
        transfer_vars(self,newM)
        return newM

    def remove(self, component):
        """Remove items by looking up column_data/row_data"""
        if type(component) is list:
            for item in component:
                if item in self.column_data: 
                    delcol = self.column_data.index(item)
                    newM = np.delete(self, delcol, axis=1).view(type(self))
                    transfer_vars(self,newM)
                    newM.column_data = list(np.delete(self.column_data, delcol))
                else:
                    delrow = self.row_data.index(item)
                    newM = np.delete(self, delrow, axis=0).view(type(self))
                    transfer_vars(self, newM)
                    # newM.column_data = self.column_data
                    # newM.model = self.model
                    newM.row_data = list(np.delete(self.row_data, delrow))

        else:
            item = component
            if item in self.column_data: 
                delcol = self.column_data.index(item)
                newM = np.delete(self, delcol, axis=1).view(type(self))
                newM.row_data = self.row_data
                newM.model = self.model
                newM.column_data = list(np.delete(self.column_data, delcol))
                # try: newM.rel = self.rel
                # except: pass
            else:
                delrow = self.row_data.index(item)
                newM = np.delete(self, delrow, axis=0).view(type(self))
                newM.column_data = self.column_data
                newM.model = self.model
                newM.row_data = list(np.delete(self.row_data, delrow))
                # try: newM.rel = self.rel
                # except: pass
        return newM

    def get(self, row_name, col_name):
        idxr = np.where(self.row_data == row_name)
        idxc = np.where(self.column_data == col_name)
        return self[idxr, idxc]


    def rows(self, component):
        rows = [self.row_data.index(item) for item in component]
        newM = self[rows,:]
        newM.model = self.model
        newM.column_data = self.column_data
        newM.row_data = list(np.array(self.row_data)[rows])
        return newM

    def del_zeros(self):
        """Delete rows and columns of a matrix with all zeros"""
        delrows = np.where(~self.any(axis=1))[0]
        delcols = np.where(~self.any(axis=0))[0]

        newM = np.delete(self, delcols, axis=1).view(type(self))
        newM = np.delete(newM, delrows, axis=0).view(type(self))
        transfer_vars(self, newM)

        newM.column_data = list(np.delete(self.column_data, delcols))
        newM.row_data = list(np.delete(self.row_data, delrows))
        return newM

    def add_cols(self, component):
        if "colinear" in component:
            vertical = [elem for elem in self.model.elems if elem.Dx==0.0]
            other = set({})
            for elem in self.model.elems:
                try: other.add(elem.Dy/elem.Dx)
                except ZeroDivisionError: pass
            return vertical, other
        if type(component) is list:
            #ASSUMES component IS A LIST OF COLUMN INDICES
            delcols = [self.column_data.index(item) for item in component[1:len(component)]]
            i0 = self.column_data.index(component[0])
            newM = np.delete(self, delcols, axis=1).view(type(self))
            for col in delcols:
                newM[:,i0] += self[:,col]
            newM.row_data = self.row_data
            newM.model = self.model
            newM.column_data = list(np.delete(self.column_data, delcols))
            # try: newM.rel = self.rel
            # except: pass
            return newM

    def add_rows(self, component):
        if type(component) is list:
            delrows = [self.row_data.index(item) for item in component[1:len(component)]]
            i0 = self.row_data.index(component[0])
            newA = np.delete(self, delrows, axis=0).view(type(self))
            for row in delrows:
                newA[i0,:] += self[row,:]
            newA.column_data = self.column_data
            newA.model = self.model
            newA.row_data = list(np.delete(self.row_data, delrows))
            return newA

class row_vector (Structural_Vector):
    def __new__(cls, Matrix):
        V = np.zeros((len(Matrix.row_data)))
        return np.asarray(V).view(cls)

    def __init__(self, Matrix):
        self.tag = "Y"
        self.matrix = Matrix 
        self.row_data = Matrix.row_data

class column_vector (Structural_Vector):

    def __new__(cls, Matrix, Vector=None):
        V = np.zeros((len(Matrix.column_data)))
        return np.asarray(V).view(cls)

    def __init__(self, Matrix, Vector=None):
        self.tag = "X"
        self.matrix = Matrix 
        self.row_data = Matrix.column_data
        if Vector is not None:
            for key in Vector.row_data:
                self.set_item(key, Vector.rows([key]))


class Static_matrix (Structural_Matrix):
    """B_MATRIX static matrix of structural model with 2d/3d truss and 2d frame elements
    the function forms the static matrix B for all degrees of freedom and
    all basic forces of the structural model specified in data structure MODEL;
    the function is currently limited to 2d/3d truss and 2d frame elements

    Parameters
    ---------------

    model: ema.Model object

    Partitions
    =========================================================================================

    - B.f  : nf x ntq

    - B.c  : nf x nq

    - B.fc : nf x nq

    - B.i  : ni x nq

    - B.x  : nx x nq

    where:

    - ni: number of primary (non-redundant) forces.
    - nq: number of total, continuous forces.
    - nx: number of redundant forces.

    """

    ranges = {
        'f': 'free dof rows',
        'i': 'primary force columns',
        'x': 'redundant force columns',
        'c': 'continuous force columns (removes releases)',
        'd': 'reaction force columns',
    }

    def __new__(cls, model, matrix=None, rng=None):
        fullnq = sum([len(elem.rel) for elem in model.elems])
        B = np.zeros((model.nt,fullnq))
        ci = 0
        for elem in model.elems:
            dofs = elem.dofs
            bg = elem.bg_matrix()
            nq = np.size(bg,1)
            for j,dof in enumerate(dofs):
                B[int(dof)-1,ci:ci+nq] = bg[j,:]
            ci = ci+nq
        input_array = B
        return np.asarray(input_array).view(cls)

    def __init__(self, model, matrix=None, rng=None):
        if rng is None:
            self.rng = None
        self.model = model
        self.row_data = np.array([str(dof) for dof in range(1, model.nt+1)])
        self.column_data = np.array([elem.tag+'_'+key for elem in model.elems for key in elem.rel.keys()])
        if matrix is not None:
            fullnq = sum([len(elem.rel) for elem in model.elems])
            self[:,:] = np.zeros((model.nt,fullnq))
            for idxr, rw in enumerate(self.row_data):
                if rw in matrix.row_data:
                    for idxc, cl in enumerate(self.column_data):
                        if cl in matrix.column_data:
                            self[idxr,idxc] = matrix.get(rw,cl)

    def __matmul__(self, Vector):
        if type(Vector) is nForce_vector:
            vect = np.matmul(self, Vector).view(iForce_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif type(Vector) is iForce_vector:
            vect = np.matmul(self, Vector).view(nForce_vector) 
            vect.row_data = self.row_data
            vect.matrix = self
        else:
            vect = Structural_Matrix.__matmul__(self, Vector)
        return vect

    @property
    def f(self):
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) > self.model.nf]

        newB = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self, newB)
        # newB.model = self.model
        newB.row_data = list(np.delete(self.row_data, delrows))
        # newB.column_data = self.column_data
        return newB

    @property
    def i(self):
        """Removes rows of B_matrix corresponding to primary (non-redundant) forces"""
        Bf = self.f
        idx_i = self.model.idx_i
        newB = Bf[:,idx_i]
        transfer_vars(Bf, newB)
        newB.column_data = Bf.column_data[idx_i]
        return newB

    @property
    def d(self):
        """Removes rows corresponding to free dofs"""
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) <= self.model.nf]
        newB = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self,newB)
        newB.row_data = list(np.delete(self.row_data, delrows))
        return newB

    @property
    def c(self):

        Bf = self.f
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        delcols = [Bf.column_data.index(tag) for tag in tags]
        newB = np.delete(Bf, delcols, axis=1).view(type(self))
        transfer_vars(Bf, newB)
        newB.column_data = list(np.delete(Bf.column_data, delcols))
        return newB

    @property
    def o(self):
        """Remove columns corresponding to element force releases, then delete zeros"""
        Bf = self.f
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delcols = [Bf.column_data.index(tag) for tag in tags]
        newB = np.delete(Bf, delcols, axis=1).view(type(self))
        transfer_vars(Bf, newB)
        newB.column_data = list(np.delete(Bf.column_data, delcols))
        newB = newB.del_zeros()
        return newB

    @property
    def fc(self):
        return self.f.c

    @property
    def x(self):
        """Removes rows of B_matrix corresponding to primary (non-redundant) forces"""
        idx_x = self.model.idx_x
        newB = self[:,idx_x]
        transfer_vars(self, newB)
        newB.column_data = self.column_data[idx_x]
        return newB

    @property
    def barxi(self):
        Bx = self.f.x
        Bbarxi = self.bari @ -Bx
        Bbarxi.column_data = Bx.column_data
        return Bbarxi

    @property
    def barx(self):
        nQ = len(self.column_data)
        nx = len(self.model.redundants)

        Bbarxi = self.barxi

        Bbarx = Structural_Matrix(np.zeros((nQ,nx)))
        transfer_vars(self, Bbarx)
        Bbarx.column_data = Bbarxi.column_data
        Bbarx.row_data = self.column_data

        for idxc, cl in enumerate(Bbarx.column_data):
            for idxr, rw in enumerate(Bbarx.row_data):
                if rw in Bbarxi.row_data:
                    Bbarx[idxr,idxc] = Bbarxi.get(rw,cl)
                elif cl==rw:
                        Bbarx[idxr, idxc] = 1.

        return Bbarx

    @property
    def bari(self):
        return self.i.del_zeros().inv

    @property
    def c0(self):
        Bf = self.f
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        delcols = [Bf.column_data.index(tag) for tag in tags]
        newB = Bf 
        for col in delcols:
            newB[:,col] = [0.]*len(Bf[:,0])
        transfer_vars(Bf, newB)
        newB.column_data = list(np.delete(Bf.column_data, delcols))
        return newB

class nStatic_matrix(Structural_Matrix):
    ranges = {
        'f': 'free dof rows',
        'i': 'primary force columns',
        'x': 'redundant force columns',
        'c': 'continuous force columns (removes releases)',
        'd': 'reaction force columns',}
    def __new__(cls, arry, model, rcdata):
        return np.asarray(arry).view(cls)

    def __init__(self, arry, model, rcdata):
        self.model = model
        self.row_data = rcdata[0]
        self.column_data = rcdata[1]

    def __matmul__(self, Vector):
        if type(Vector) is nForce_vector:
            vect = np.matmul(self, Vector).view(iForce_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif isinstance(Vector, iForce_vector):
            vect = np.matmul(self, Vector).view(nForce_vector) 
            vect.row_data = self.row_data
            vect.matrix = self
        else:
            vect = Structural_Matrix.__matmul__(self, Vector)
        return vect

    @property
    def f(self):
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) > self.model.nf]

        newB = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self, newB)
        newB.row_data = np.delete(self.row_data, delrows)
        return newB

    @property
    def i(self):
        """Removes rows of B_matrix corresponding to redundant forces"""
        # reducedB = self.f.c.del_zeros()
        reducedB = self.f.c
        rdts = reducedB.model.redundants
        tags = [q.elem.tag + '_'+str(q.nature) for q in rdts]
        delcols = [reducedB.column_data.index(tag) for tag in tags]
        newB = np.delete(reducedB, delcols, axis=1).view(Static_matrix)
        transfer_vars(reducedB, newB)
        newB.column_data = list(np.delete(reducedB.column_data, delcols))
        return newB

    @property
    def d(self):
        """Removes rows corresponding to free dofs"""
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) <= self.model.nf]
        newB = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self,newB)
        newB.row_data = list(np.delete(self.row_data, delrows))
        return newB

    @property
    def c(self):
        """Removes columns corresponding to element hinges/releases"""
        Af = self.f
        idx_c = self.model.idx_c
        newA = Af[:,idx_c]
        transfer_vars(Af, newA)
        newA.column_data = Af.column_data[idx_c]
        return newA

    @property
    def o(self):
        """Remove columns corresponding to element force releases, then delete zeros"""
        Bf = self.f
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delcols = [Bf.column_data.index(tag) for tag in tags]
        newB = np.delete(Bf, delcols, axis=1).view(type(self))
        transfer_vars(Bf, newB)
        newB.column_data = list(np.delete(Bf.column_data, delcols))
        newB = newB.del_zeros()
        return newB

    @property
    def fc(self):
        return self.f.c

    @property
    def x(self):
        """Removes rows of B_matrix corresponding to primary (non-redundant) forces

        """
        rdts = self.model.redundants
        tags = [q.elem.tag + '_'+str(q.nature) for q in rdts]
        cols = [self.column_data.index(tag) for tag in tags]
        newB = self[:,cols]
        transfer_vars(self, newB)
        newB.column_data = np.array([self.column_data[col] for col in cols])
        # newB.row_data = self.row_data
        # newB.model = self.model
        return newB

    @property
    def barxi(self):
        Bx = self.f.x
        Bbarxi = self.bari @ -Bx
        Bbarxi.column_data = Bx.column_data
        return Bbarxi

    @property
    def barx(self):
        nQ = len(self.column_data)
        nx = len(self.model.redundants)

        Bbarxi = self.barxi

        Bbarx = Structural_Matrix(np.zeros((nQ,nx)))
        transfer_vars(self, Bbarx)
        Bbarx.column_data = Bbarxi.column_data
        Bbarx.row_data = self.column_data
        for idxc, cl in enumerate(Bbarx.column_data):
            for idxr, rw in enumerate(Bbarx.row_data):
                if rw in Bbarxi.row_data:
                    Bbarx[idxr,idxc] = Bbarxi.get(rw,cl)
                elif cl==rw:
                        Bbarx[idxr, idxc] = 1.

        return Bbarx

    @property
    def bari(self):
        return self.i.del_zeros().inv

    @property
    def c0(self):
        Bf = self.f
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delcols = [Bf.column_data.index(tag) for tag in tags]
        newB = Bf 
        for col in delcols:
            newB[:,col] = [0.]*len(Bf[:,0])
        transfer_vars(Bf, newB)
        newB.column_data = list(np.delete(Bf.column_data, delcols))
        # newB.row_data = self.row_data
        # newB.rel = self.rel
        return newB

def B_matrix(model, matrix=None, rng=None):
    """Returns a Static_matrix object"""
    return Static_matrix(model, matrix, rng)

def nB_matrix(model):
    """Returns a Static_matrix object"""
    fullnq = sum([len(elem.rel) for elem in model.elems])
    B = np.zeros((model.nt, fullnq))
    ci = 0
    for elem in model.elems:
        eid = elem.dofs
        bg = elem.bg_matrix()
        nq = np.size(bg,1)
        for j,eidi in enumerate(eid):
            B[int(eidi)-1,ci:ci+nq] = bg[j,:]
        ci = ci+nq
    input_array = B
    # return np.asarray(input_array).view(nStatic_matrix)
    matrix =  np.asarray(input_array)
    return nStatic_matrix(model, matrix, rng)

def Bh_matrix(model):
    """Returns a Static_matrix object"""
    fullnq = sum([len(elem.rel) for elem in model.elems])
    B = np.zeros((model.nt, fullnq))
    ci = 0
    for elem in model.elems:
        eid = elem.dofs
        bg = elem.bg_matrix(Roption=True)
        nq = np.size(bg,1)
        for j,eidi in enumerate(eid):
            B[int(eidi)-1,ci:ci+nq] = bg[j,:]
        ci = ci+nq

    row_data = np.array([str(dof) for dof in range(1, model.nt+1)])
    column_data = np.array([elem.tag+'_'+key
                            for elem in model.elems
                            for key in elem.rel.keys()])
    rcdata = (row_data, column_data)
    return nStatic_matrix(B, model, rcdata)

class Kinematic_matrix(Structural_Matrix):
    """Class for the kinematic matrix of a structural model with 2d/3d truss and 2d frame elements
    the function forms the kinematic matrix A for all degrees of freedom and
    all element deformations of the structural model specified in data structure MODEL 
    the function is currently limited to 2d/3d truss and 2d frame elements

    Returns
    ---------

    Kinematic matrix

    """

    ranges = {
        'f': 'free dof columns',
        'i': 'primary force/deformation rows',
        'x': 'redundant force/deformation rows',
        'd': 'reaction force/deformation rows',
        'c': 'continuous (hingeless) force/deformation rows'
    }

    def __new__(cls, model, matrix=None,rng=None):
        A  = np.zeros((sum(model.nv),model.nt))
        ri = 0
        for elem in model.elems:
            eid = elem.dofs
            ag = elem.ag()
            nv = len(elem.v)
            for j, eidi in enumerate(eid):
                A[ri:ri+nv, int(eidi)-1] = ag[:,j]
            ri = ri+nv
        input_array = A
        return np.asarray(input_array).view(cls)

    def __init__(self, model, matrix=None,rng=None):
        if rng is None:
            self.rng = None
        self.model = model
        self.column_data = np.array([str(dof) for dof in range(1,model.nt+1)])
        self.row_data = np.array([elem.tag+'_'+key for elem in model.elems for key in elem.v.keys()])
        self.basic_deformations = np.array([v for elem in model.elems for v in elem.basic_deformations])

        self.idx_h = []
        if matrix is not None:
            fullnq = sum([len(elem.rel) for elem in model.elems])
            self[:,:] = np.zeros((model.nt,fullnq))
            for idxr, rw in enumerate(self.row_data):
                if rw in matrix.row_data:
                    for idxc, cl in enumerate(self.column_data):
                        if cl in matrix.column_data:
                            self[idxr,idxc] = matrix.get(rw,cl)

    def __matmul__(self, Vector):
        if isinstance(Vector, Deformation_vector):
            # print('a')
            vect = np.matmul(self, Vector).view(Displacement_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif isinstance(Vector, Displacement_vector):
            vect = np.matmul(self, Vector).view(Deformation_vector) 
            vect.row_data = self.row_data
            vect.matrix = self
        else:
            vect = Structural_Matrix.__matmul__(self, Vector)

        return vect


    def combine(self, component):
        if "colinear" in component:
            vertical = [elem for elem in self.model.elems if elem.Dx==0.0]
            other = set({})
            for elem in self.model.elems:
                try: other.add(elem.Dy/elem.Dx)
                except ZeroDivisionError: pass
            return vertical, other

        if type(component) is list:
            ## TO BE DEPRECATED
            #ASSUMES component IS A LIST OF COLUMN INDICES
            delcols = [self.column_data.index(item) for item in component[1:len(component)]]
            i0 = self.column_data.index(component[0])
            newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
            for col in delcols:
                newA[:,i0] += self[:,col]
            newA.row_data = self.row_data
            newA.model = self.model
            newA.column_data = list(np.delete(self.column_data, delcols))
            # newA.rel = self.rel
            return newA

    @property
    def f(self):
        """Removes columns corresponding to fixed dofs"""
        delcols = [idx for idx, dof in enumerate(self.column_data) if int(dof) > self.model.nf]
        newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
        transfer_vars(self, newA)
        newA.column_data = list(np.delete(self.column_data, delcols))
        return newA

    @property
    def i(self):
        """Removes rows corresponding to redundant forces"""
        Afc = self.f.c
        rdts = self.model.redundants
        tags = [q.elem.tag +'_'+ str(q.nature) for q in rdts]
        delrows = [Afc.row_data.index(tag) for tag in tags]
        newA = np.delete(Afc, delrows, axis=0).view(Kinematic_matrix)
        transfer_vars(Afc,newA)
        newA.row_data = list(np.delete(Afc.row_data, delrows))
        return newA

    @property
    def d(self):
        """Removes columns corresponding to free dofs"""
        delcols = [idx for idx, dof in enumerate(self.column_data) if int(dof) <= self.model.nf]
        newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
        transfer_vars(self,newA)
        newA.column_data = list(np.delete(self.column_data, delcols))
        return newA


    @property
    def c(self):
        """Removes rows corresponding to element hinges/releases"""
        Af = self.f
        idx_c = self.model.idx_c
        newA = Af[idx_c,:]
        transfer_vars(Af, newA)
        newA.row_data = Af.row_data[idx_c]
        return newA

    @property
    def c0(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = np.where(np.isin(Af.row_data,tags))
        # delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA

    @property
    def o(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        newA = newA.del_zeros()
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA

    @property
    def e(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        newA = newA.del_zeros()
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA


class nKinematic_matrix (Structural_Matrix):
    """Class for the kinematic matrix of a structural model with 2d/3d truss and 2d frame elements
    the function forms the kinematic matrix A for all degrees of freedom and
    all element deformations of the structural model specified in data structure MODEL 
    the function is currently limited to 2d/3d truss and 2d frame elements

    Returns
    ---------
    Kinematic matrix

    """

    ranges = {
        'f': 'free dof columns',
        'i': 'primary force/deformation rows',
        'x': 'redundant force/deformation rows',
        'd': 'reaction force/deformation rows',
        'c': 'continuous (hingeless) force/deformation rows'
    }

    def __new__(cls, model, matrix=None,rng=None):
        A  = np.zeros((sum(model.nv),model.nt))
        ri = 0
        for elem in model.elems:
            eid = elem.dofs
            ag = elem.ag()
            nv = len(elem.v)
            for j, eidi in enumerate(eid):
                A[ri:ri+nv, int(eidi)-1] = ag[:,j]
            ri = ri+nv
        input_array = A
        return np.asarray(input_array).view(cls)

    def __new__(cls, arry, model, rcdata):
        return np.asarray(arry).view(cls)

    def __init__(self, arry, model, rcdata):
        self.model = model
        self.row_data = rcdata[0]
        self.column_data = rcdata[1]
        self.basic_deformations = np.array([v for elem in model.elems for v in elem.basic_deformations])

        self.idx_h = []


    def __matmul__(self, Vector):
        if isinstance(Vector, Deformation_vector):
            # print('a')
            vect = np.matmul(self, Vector).view(Displacement_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif isinstance(Vector, Displacement_vector):
            # print('b')
            vect = np.matmul(self, Vector).view(Deformation_vector)
            # vect = np.matmul(self, Vector).view(Structural_Vector)  
            vect.row_data = self.row_data
            vect.matrix = self
        else:
            # print('else')
            vect = Structural_Matrix.__matmul__(self, Vector)
        return vect

    def combine(self, component):
        if "colinear" in component:
            vertical = [elem for elem in self.model.elems if elem.Dx==0.0]
            other = set({})
            for elem in self.model.elems:
                try: other.add(elem.Dy/elem.Dx)
                except ZeroDivisionError: pass
            return vertical, other

        if type(component) is list:
            ## TO BE DEPRECATED
            #ASSUMES component IS A LIST OF COLUMN INDICES
            delcols = [self.column_data.index(item) for item in component[1:len(component)]]
            i0 = self.column_data.index(component[0])
            newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
            for col in delcols:
                newA[:,i0] += self[:,col]
            newA.row_data = self.row_data
            newA.model = self.model
            newA.column_data = list(np.delete(self.column_data, delcols))
            # newA.rel = self.rel
            return newA

    @property
    def f(self):
        """Removes columns corresponding to fixed dofs"""
        delcols = [idx for idx, dof in enumerate(self.column_data) if int(dof) > self.model.nf]
        newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
        transfer_vars(self, newA)
        newA.column_data = list(np.delete(self.column_data, delcols))
        return newA

    @property
    def i(self):
        """Removes rows corresponding to redundant forces"""
        Afc = self.f.c
        rdts = self.model.redundants
        tags = [q.elem.tag +'_'+ str(q.nature) for q in rdts]
        delrows = [Afc.row_data.index(tag) for tag in tags]
        newA = np.delete(Afc, delrows, axis=0).view(Kinematic_matrix)
        transfer_vars(Afc,newA)
        newA.row_data = list(np.delete(Afc.row_data, delrows))
        return newA

    @property
    def d(self):
        """Removes columns corresponding to free dofs"""
        delcols = [idx for idx, dof in enumerate(self.column_data) if int(dof) <= self.model.nf]
        newA = np.delete(self, delcols, axis=1).view(Kinematic_matrix)
        transfer_vars(self,newA)
        newA.column_data = list(np.delete(self.column_data, delcols))
        return newA


    @property
    def c(self):
        """Removes rows corresponding to element hinges/releases"""
        Af = self.f
        idx_c = self.model.idx_c
        newA = Af[idx_c,:]
        transfer_vars(Af, newA)
        newA.row_data = Af.row_data[idx_c]
        return newA

    @property
    def c0(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA

    @property
    def o(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        newA = newA.del_zeros()
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA

    @property
    def e(self):
        Af = self.f
        n_col = len(Af.T)
        tags = [elem.tag + "_" + rel for elem in self.model.elems for rel in elem.rel if elem.rel[rel]]
        # delcols = [idx for idx, rel in enumerate(self.model.rel) if rel==1]
        delrows = [Af.row_data.index(tag) for tag in tags]
        newA = Af 
        for rw in delrows:
            newA[rw,:] = [0.]*n_col
        newA = newA.del_zeros()
        transfer_vars(Af, newA)
        # newA.row_data = list(np.delete(Af.row_data, delrows))
        return newA

def A_matrix(Domain, matrix=None):
    """Returns a Kinematic_matrix object"""
    return Kinematic_matrix(Domain,matrix)


class Diag_matrix(Structural_Matrix):
    """Block diagonal matrix of element flexibility/stiffness matrices for structural model


    this class represents the block diagonal matrix of element flexibility or stiffness matrices
    for a structural model.

    """
    tag = 'F'
    def __new__(cls, arry, rc_data, model):
        basic_forces = rc_data
        arry = np.asarray(arry).view(cls)
        arry.basic_forces = basic_forces
        arry.model = model
        # arry.column_data = [elem.tag+'_'+key for elem in model.elems for key in elem.rel.keys()]
        arry.column_data = arry.row_data = np.array([q.elem.tag+'_'+str(q.number) for q in basic_forces])
        # arry.row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.rel.keys()]
        return arry

    def __init__(self, arry, rc_data, model):
        """Parameters
        =========================================================================================
        model

        """
        self.model = model

    def __matmul__(self, Vector):
        if type(Vector) is Deformation_vector:
            vect = np.matmul(self, Vector).view(iForce_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif isinstance(Vector, iForce_vector):
            vect = np.matmul(self, Vector).view(Deformation_vector)
            # vect.part = 'continuous'
            vect.row_data = self.row_data
            vect.matrix = self

        else:
            vect = Structural_Matrix.__matmul__(self, Vector)
        return vect

    @property
    def c(self):
        """Removes columns corresponding to element hinges/releases"""
        idx_c = self.model.idx_c
        newA = self[:,idx_c]
        newA = newA[idx_c,:]
        transfer_vars(self, newA)
        newA.basic_forces = self.basic_forces[idx_c]
        newA.column_data = self.column_data[idx_c]
        newA.row_data = self.row_data[idx_c]
        return newA

def Fs_matrix(model, Roption=True):
    """Returns a Flexibility_matrix object"""  
    if Roption:
        f  = np.array([elem.f_matrix(Roption) for elem in model.elems])
        basic_forces = np.array([q for elem in model.elems for q in elem.basic_forces if not q.rel]) 
    else:
        f  = np.array([elem.f_matrix(Roption) for elem in model.elems])
        basic_forces = np.array([q for elem in model.elems for q in elem.basic_forces]) 
    Fs = scipy.linalg.block_diag(*f) 

    basic_forces = basic_forces
    Fs = Diag_matrix(Fs, basic_forces, model)
    return Fs

def Ks_matrix(model):
    """Returns a Flexibility_matrix object"""  
    k  = np.array([elem.k_matrix() for elem in model.elems])
    Ks = scipy.linalg.block_diag(*k)  
    basic_forces = np.array([q for elem in model.elems for q in elem.basic_forces])
    Ks = Diag_matrix(Ks, basic_forces, model)
    return Ks


class Stiffness_matrix (Structural_Matrix):
    """...
    Parameters
    =========================================================================================
    model

    -----------------------------------------------------------------------------------------
    """
    tag = 'K'
    def __new__(cls, arry, model, Roption=None):

        input_array = np.asarray(arry).view(cls)
        input_array.model = model
        input_array.column_data = [str(dof) for dof in range(1, model.nt+1)]
        input_array.row_data = ['P_{'+str(dof)+'}' for dof in range(1, model.nt+1)]

        return input_array

    def __init__(self, arry, model, Roption=None):
        self.subs = [None]
        pass


    def __matmul__(self, Vector):
        if type(Vector) is Displacement_vector:
            vect = np.matmul(self, Vector).view(nForce_vector)
            vect.row_data = self.row_data
            vect.matrix = self
        elif type(Vector) is nForce_vector:
            vect = np.matmul(self, Vector).view(Displacement_vector)
            vect.row_data = self.row_data
            vect.model = self.model
        else:
            vect = Structural_Matrix.__matmul__(self, Vector)
        return vect

    @property
    def f(self):
        delrows = [idx for idx, dof in enumerate([str(dof) for dof in range(1, self.model.nt+1)]) if int(dof) > self.model.nf]

        newK = np.delete(self, delrows, axis=0)
        newK = np.delete(newK, delrows, axis=1).view(type(self))
        transfer_vars(self, newK)

        newK.row_data = list(np.delete(self.row_data, delrows))
        newK.column_data = list(np.delete(self.column_data, delrows))
        return newK

def K_matrix(Model):
    """Returns a Stiffness_matrix object"""
    K = np.zeros((Model.nt, Model.nt))
    for elem in Model.elems:
        ke = elem.ke_matrix()
        for i, dof in enumerate(elem.dofs):
            for j, doff in enumerate(elem.dofs):
                K[int(dof)-1, int(doff)-1] += ke[i, j]

    return Stiffness_matrix(K, Model, Roption=None)



def Kt_matrix(Model, State):
    """Returns a Stiffness_matrix object"""
    K = np.zeros((Model.nt, Model.nt))
    for elem in Model.elems:
        kt = elem.kt_matrix(State)
        for i, dof in enumerate(elem.dofs):
            for j, doff in enumerate(elem.dofs):
                K[int(dof)-1, int(doff)-1] += kt[i, j]

    return Stiffness_matrix(K, Model, Roption=None)



class Displacement_vector(column_vector):
    tag = 'U'
    def __new__(cls, Kinematic_matrix, Vector=None):
        U = np.zeros((len(Kinematic_matrix.column_data)))
        return np.asarray(U).view(cls)

    def __init__(self, Kinematic_matrix, Vector=None):
        self.matrix = Kinematic_matrix 
        self.row_data = Kinematic_matrix.column_data
        self.subs = [None]
        if Vector is not None:
            for key in Vector.row_data:
                if key in self.row_data:
                    self.set_item(key, Vector.rows([key]))

    @property
    def f(self):
        """Removes rows corresponding to fixed dofs"""
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) > self.model.nf]
        newU = np.delete(self, delrows, axis=0).view(Displacement_vector)
        newU.row_data = list(np.delete(self.row_data, delrows))
        # newU.matrix = self.matrix
        newU.model = self.model
        return newU

class nDisplacement_vector(Structural_Vector):
    tag = 'U'
    def __new__(cls, arry, model, row_data, Vector=None):
        input_array = np.asarray(arry).view(cls)
        return input_array

    def __init__(self, arry, model, row_data, Vector=None):
        self.subs = [None]
        self.model = model
        self.row_data = row_data

        if Vector is not None:
            for key in Vector.row_data:
                if key in self.row_data:
                    self.set_item(key, Vector.rows([key]))

    @property
    def f(self):
        """Removes rows corresponding to fixed dofs"""
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) > self.model.nf]
        newU = np.delete(self, delrows, axis=0).view(nDisplacement_vector)
        newU.row_data = list(np.delete(self.row_data, delrows))
        newU.model = self.model
        return newU

def U_vector(model, vector=None):
    """Returns a Displacement_vector object"""   
    U = np.zeros(model.nt)
    row_data = [str(dof) for dof in range(1,model.nt+1)]
    U = nDisplacement_vector(U, model, row_data)

    if vector is not None:
        if len(vector)==len(U):
            U[:,0] = vector[:]
        else:
            for key in vector.row_data:
                if key in U.row_data:
                    U.set_item(key, vector.rows([key]))
    return U


class iForce_vector(Structural_Vector):
    tag = 'Q'
    def __new__(cls, arry, model, row_data, Vector=None):
        return np.asarray(arry).view(cls)

    def __init__(self, arry, model, row_data, Vector=None):
        self.model = model
        self.subs = [None]
        self.row_data = row_data
        if Vector is not None:
            for key in Vector.row_data:
                if key in self.row_data:
                    self.set_item(key, Vector.get(key))
  #<
    @property
    def i(self):
        """Removes rows corresponding to redundant forces"""
        rdts = self.model.redundants

        tags = [q.elem.tag + str(q.nature) for q in rdts]
        delrows = [self.row_data.index(tag) for tag in tags]
        newQ = np.delete(self, delrows, axis=0).view(iForce_Vector)
        transfer_vars(self, newQ)
        newQ.subs.append('primary')
        newQ.row_data = list(np.delete(self.row_data, delrows))
        return newQ

    @property
    def c(self):
        """Remove rows corresponding to element hinges/releases"""
        idx_c = self.model.idx_c
        newQ = self[idx_c]
        transfer_vars(self, newQ)
        newQ.row_data = self.row_data[idx_c]
        return newQ

    @property
    def x(self):
        """Remove rows of corresponding to primary forces"""
        rdts = self.model.redundants
        tags = [q.elem.tag + '_'+str(q.nature) for q in rdts]
        rows = [self.row_data.index(tag) for tag in tags]
        newV = self[rows]
        newV.row_data = [self.row_data[row] for row in rows]
        transfer_vars(self, newV)
        return newV

def Q_vector(model, vector=None):
    """Returns a iForce_vector object"""   

    arry = np.zeros((model.nQ,1))
    row_data = np.array([elem.tag+'_'+key for elem in model.elems for key in elem.rel.keys()])
    Q = iForce_vector(arry, model, row_data)
    if vector is not None:
        if len(vector)==len(Q):
            Q[:,0] = vector[:]
        else:
            for key in vector.row_data:
                if key in Q.row_data:
                    Q.set_item(key, vector.rows([key]))
    return Q

def Q0_vector(model):
    """Returns a vector of initial element forces"""   
    arry = np.concatenate([elem.q0_vector() for elem in model.elems])
    row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.q.keys()] 
    return iForce_vector(arry, model, row_data)

def Qpl_vector(model):
    """Returns a vector of element plastic capacities"""

    Qp_pos = [elem.Qp['+'][key] for elem in model.elems for key in elem.Qp['+']]
    Qp_neg = [elem.Qp['-'][key] for elem in model.elems for key in elem.Qp['-']]
    row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.Qp['-']]

    # del_idx = np.where(~Bf.any(axis=0))[0]
    # Qp_pos = np.delete(Qp_pos, del_idx)
    # Qp_neg = np.delete(Qp_neg, del_idx)
    # row_data = np.delete(row_data, del_idx)
    column_data = ['Q_{pl}^+', 'Q_{pl}^-']

    Qpl = nKinematic_matrix(np.array([Qp_pos, Qp_neg]).T, model, (row_data, column_data))
    return Qpl

def Qp_vector(model):
    """Returns a vector of element plastic capacities"""
    B = B_matrix(model)
    Bf = B.f
    Qp_pos   = [elem.Qp['+'][key] for elem in model.elems for key in elem.Qp['+']]
    Qp_neg   = [elem.Qp['-'][key] for elem in model.elems for key in elem.Qp['-']]
    row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.Qp['-']]

    del_idx = np.where(~Bf.any(axis=0))[0]
    Qp_pos = np.delete(Qp_pos, del_idx)
    Qp_neg = np.delete(Qp_neg, del_idx)
    row_data = np.delete(row_data, del_idx)
    column_data = ['Q_{pl}^+', 'Q_{pl}^-']

    Qpl = nStatic_matrix(np.array([Qp_pos, Qp_neg]).T, model, (row_data, column_data))
    return Qpl

class nForce_vector(Structural_Vector):
    tag = 'P'
    def __new__(cls, arry, model, row_data, Vector=None):
        return np.asarray(arry).view(cls)

    def __init__(self, arry, model, row_data, Vector=None):
        self.model = model
        self.subs = [None]
        self.row_data = row_data
        if Vector is not None:
            for key in Vector.row_data:
                if key in self.row_data:
                    self.set_item(key, Vector.get(key))


    @property
    def f(self):
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) > self.model.nf]

        newP = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self, newP)
        # newB.model = self.model
        newP.row_data = list(np.delete(self.row_data, delrows))
        # newB.column_data = self.column_data
        return newP

    @property
    def d(self):
        """Removes rows corresponding to free dofs"""
        delrows = [idx for idx, dof in enumerate(self.row_data) if int(dof) <= self.model.nf]
        newP = np.delete(self, delrows, axis=0).view(type(self))
        newP.model = self.model
        newP.row_data = list(np.delete(self.row_data, delrows))
        newP.subs.append('reactions')
        return newP

  #<

def P_vector(model, vector=None):
    P = np.zeros(model.nt)
    for node in model.nodes:
        p = node.p_vector()
        for i, dof in enumerate(node.dofs):
            P[int(dof)-1] += p[i]
    row_data = [str(dof) for dof in range(1, model.nt+1)]
    return nForce_vector(P, model, row_data)



# def P0_vector(model):
#     """Returns a _ object"""   
#     arry = np.concatenate([elem.p0_vector() for elem in model.elems]) 
#     row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.v.keys()] 
#     return nForce_vector(arry, model, row_data)

def P0_vector(model):
    P = np.zeros(model.nt)
    for elem in model.elems:
        dofs = elem.dofs
        if hasattr(elem, 'q0_vector'):
            p0 = elem.bg_matrix()@elem.q0_vector()
        else:
            p0 = np.zeros(len(dofs))

        for i,df in enumerate(dofs):
            P[int(df)-1] +=  p0[i]

    row_data = [str(dof) for dof in range(1, model.nt+1)]
    return nForce_vector(P, model, row_data)

def Pw_vector(model):
    P = np.zeros(model.nt)
    for elem in model.elems:
        dofs = elem.dofs
        if len(dofs)==6:
            P[int(dofs[0])-1] +=  elem.w['y']*elem.L/2*elem.sn
            P[int(dofs[1])-1] += -elem.w['y']*elem.L/2*elem.cs
            P[int(dofs[3])-1] +=  elem.w['y']*elem.L/2*elem.sn
            P[int(dofs[4])-1] += -elem.w['y']*elem.L/2*elem.cs
        else:
            pw = elem.pw_vector()
            for i,df in enumerate(dofs):
                P[int(df)-1] +=  pw[i]


    row_data = [str(dof) for dof in range(1, model.nt+1)]
    return nForce_vector(P, model, row_data)


class Deformation_vector(Structural_Vector):
    tag = "V"
    def __new__(cls, arry, model, row_data, Vector=None):
        input_array = np.asarray(arry).view(cls)
        return input_array

    def __init__(self, arry, model, row_data, Vector=None):
        self.model = model
        self.row_data = row_data

    @property
    def c(self):
        """Removes rows corresponding to element hinges/releases"""
        idx_c = self.model.idx_c
        newQ = self[idx_c]
        transfer_vars(self, newQ)
        newQ.row_data = self.row_data[idx_c]
        return newQ

    @property
    def i(self):
        """Removes rows corresponding to redundant forces"""
        rdts = self.model.redundants
        tags = [q.elem.tag + '_'+str(q.nature) for q in rdts]
        delrows = [self.row_data.index(tag) for tag in tags]
        newV = np.delete(self, delrows, axis=0).view(type(self))
        transfer_vars(self,newV)
        newV.row_data = list(np.delete(self.row_data, delrows))
        newV.subs.append('primary')
        return newV

    @property
    def x(self):
        """Removes rows of corresponding to primary forces

        """
        rdts = self.model.redundants
        tags = [q.elem.tag + '_'+str(q.nature) for q in rdts]
        rows = [self.row_data.index(tag) for tag in tags]
        newV = self[rows]
        newV.row_data = [self.row_data[row] for row in rows]
        transfer_vars(self, newV)
        return newV

def V_vector(model, vector=None):
    """Returns a Deformation_vector object"""
    arry = np.zeros((model.nQ,1))
    row_data = np.array([elem.tag+'_'+key for elem in model.elems for key in elem.rel.keys()])
    V = Deformation_vector(arry, model, row_data)
    if vector is not None:
        if len(vector)==len(V):
            V[:,0] = vector[:]
        else:
            for key in vector.row_data:
                if key in V.row_data:
                    V.set_item(key, vector.rows([key]))
    return V

def V0_vector(model):
    """Returns a Deformation_vector object"""
    arry = np.concatenate([elem.v0_vector() for elem in model.elems]) 
    row_data = [elem.tag+'_'+key for elem in model.elems for key in elem.v.keys()] 
    return Deformation_vector(arry, model, row_data)


def Aub_matrix(model, alpha):
    """Return the interaction upperbound matrix"""
    aub  = np.array([elem.aub_matrix(alpha) for elem in model.elems])
    A = scipy.linalg.block_diag(*aub)
    return A

