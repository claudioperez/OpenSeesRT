import numpy as np
from abc import abstractmethod

from numpy.polynomial import Polynomial
from scipy.integrate import quad
import scipy.integrate

from .matrices import Structural_Matrix, Structural_Vector


class IntForce:
    def __init__(self, elem, number):
        self.number = number
        self.elem = elem
        self.rel = False
        self.redundant = False
        self.plastic_event = None

    @property
    def tag(self):
        return self.elem.tag + '_' + str(self.number)

    def __str__(self):
        return self.tag

    def __repr__(self):
        return self.tag


class BasicLink():
    """Class implementing general geometric element methods"""

    def __init__(self, ndf, ndm, nodes):
        self.nodes = nodes
        self.ndf: int = ndf
        self.ndm: int = ndm
        self.nen = len(nodes)

    @property
    def dofs(self):
        eid = np.array([])
        for node in self.nodes:
            eid = np.append(eid, node.dofs[0:self.ndf])
        return eid

    @property
    def L(self):
        xyzi = self.nodes[0].xyz
        xyzj = self.nodes[1].xyz
        L = np.linalg.norm(xyzi-xyzj)
        return L

    @property
    def L0(self):
        xyzi = self.nodes[0].xyz0
        xyzj = self.nodes[1].xyz0
        L = np.linalg.norm(xyzi-xyzj)
        return L

    @property
    def Li(self):
        n1 = self.nodes[0]
        n2 = self.nodes[1]
        xyzi = np.array([n1.xi, n1.yi, n1.zi])
        xyzj = np.array([n2.xi, n2.yi, n2.zi])
        L = np.linalg.norm(xyzi-xyzj)
        return L


    @property
    def Dx(self):
        return self.nodes[1].x - self.nodes[0].x

    @property
    def Dy(self):
        return self.nodes[1].y - self.nodes[0].y

    @property
    def Dz(self):
        return self.nodes[1].z - self.nodes[0].z

    @property
    def sn(self):
        L = self.L
        sn = (self.nodes[1].y - self.nodes[0].y)/L
        return sn

    @property
    def cs(self):
        L = self.L
        cs = (self.nodes[1].x - self.nodes[0].x)/L
        return cs     

    @property
    def cz(self):
        L = self.L
        cz = (self.nodes[1].z - self.nodes[0].z)/L
        return cz

    def Rx_matrix(self):
        """Rotation about x

        https://en.wikipedia.org/wiki/Rotation_matrix
        """
        cs = self.cs 
        sn = self.sn
        Rx = np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0,  cs, -sn, 0.0, 0.0, 0.0],
            [0.0,  sn,  cs, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,  cs, -sn],
            [0.0, 0.0, 0.0, 0.0,  sn,  cs],
        ])

        return 0

    def Ry_matrix(self):
        "Rotation about z"
        cs = self.cs 
        sn = self.sn
        Ry = np.array([
            [ cs, 0.0,  sn],
            [0.0, 1.0, 0.0],
            [-sn, 0.0,  cs],
        ])

        return 0

    def Rz_matrix(self):
        "Rotation about z"
        cs = self.cs
        sn = self.sn
        Rz = np.array([
            [ cs, -sn, 0.0, 0.0, 0.0, 0.0],
            [ sn,  cs, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0,  cs, -sn, 0.0],
            [0.0, 0.0, 0.0,  sn,  cs, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

        return Rz


class Element(BasicLink):
    """Element parent class"""

    def __init__(self,  ndf, ndm, force_dict, nodes):
        super().__init__(  ndf, ndm, nodes)

        self.history = {}
        self.current = {}

        nq = len(force_dict)
        self.rel = {str(i): False for i in range(1, nq+1)}
        self.red = {str(i): False for i in range(1, nq+1)}
        self.q0  = {str(i): 0.  for i in range(1, nq+1)}
        self.v0  = {str(i): 0.  for i in range(1, nq+1)}
        self.e0  = {str(i): 0.  for i in range(1, nq+1)}
        self.Qp  = {'+': {str(i): 0.  for i in range(1, nq+1)}, '-':{str(i): 0.  for i in range(1, nq+1)}}

        self.basic_forces = np.array([IntForce(self, i) for i in range(1, nq+1)])
        self.basic_deformations = self.basic_forces

    @property
    def force_keys(self):
        return [self.tag+'_'+key for key in self.rel]


    def setHistoryParameter ( self, name, val ):
        self.current[name] = val


    def getHistoryParameter ( self, name ):
        return self.history[name]


    def commitHistory ( self ):
        self.history = self.current.copy()
        self.current = {}

        if hasattr( self , "mat" ):
            self.mat.commitHistory()

    def v0_vector(self):
        return np.array([0]*self.ndf*self.nn)

    def pw_vector(self):
        if all(self.w.values())==0.0:
            return np.array([0.0]*self.ndf*self.nn)

class PolyRod(Element):
    nv  = 1
    nn  = 2
    ndm = 1
    ndf = 1 # number of dofs at each node
    force_dict = {'1':0}

    def __init__(self,tag, nodes, E, A):
        super().__init__(self.ndf, self.ndm, self.force_dict, nodes)
        self.tag = tag
        self.E = E 
        self.A = A
        self.q = {'1':0}
        self.v = {'1':0}
        self.w = {'1':0.0}

        if isinstance(self.E,float):
            self.E = Polynomial([self.E])

        if isinstance(self.A,float):
            self.A = Polynomial([self.A])

    def N(self):
        L = self.L
        N1 = Polynomial([1,-1/L])
        N2 = Polynomial([0, 1/L])
        return np.array([[N1],[N2]])

    def B(self):
        N = self.N()
        return np.array([[Polynomial.deriv(Ni,1) for Ni in row] for row in N])

    def k_matrix(self):
        E = self.E
        A = self.A

        L = self.L
        B = self.B()
        k = Structural_Matrix([
            [quad(E*A*(B@B.T)[i,j],0,L)[0] for j in range(2)] for i in range(2)
            ])

        # Metadata
        k.tag = 'k'
        k.row_data = k.column_data = ['u_'+str(int(dof)) for dof in self.dofs]
        return k

    def ke_matrix(self):
        return self.k_matrix()

    def pw_vector(self):
        L = self.L
        pw = self.w['1']
        N = self.N()
        if isinstance(pw,np.polynomial.Polynomial) or isinstance(pw,float):
            # p1 = -quad(N[0]*pw,0,L)[0]
            # p2 = -quad(N[1]*pw,0,L)[0]
            p = np.array([[-quad(Ni*pw,0,L)[0] for Ni in row] for row in N])
        else:
            print('Unsupported load vector pw')

        return p

    ## Post Processing
    def localize(self,U_vector):
        dofs = [int(dof) for dof in self.dofs]
        return np.array([U_vector.get(str(dof)) for dof in dofs])

    def displx(self,U_vector):
        u = self.localize(U_vector)
        N = self.N()
        return N[0,0]*u[0] + N[1,0]*u[1]

    def strainx(self,U_vector):
        dofs = [int(dof) for dof in self.dofs]
        u = np.array([U_vector.get(str(dof)) for dof in dofs])
        B = self.B()
        return B[0,0]*u[0] + B[1,0]*u[1]

    def iforcex(self,U_vector):
        """Resisting force vector"""
        dofs = [int(dof) for dof in self.dofs]
        u = np.array([U_vector.get(str(dof)) for dof in dofs])
        B = self.B()
        P = self.E*self.A*(B[0,0]*u[0] + B[1,0]*u[1])
        return P

class Truss(Element):
    nv = 1
    nn = 2
    nq = 1
    ndm = 2
    ndf = 2 # number of dofs at each node
    force_dict = {'1':0}
    Qpl = np.zeros((2,nq))

    def __init__(self, tag, iNode, jNode, E=None, A=None, geo='lin',properties=None):
        super().__init__(self.ndf, self.ndm, self.force_dict, [iNode, jNode])
        if isinstance(properties,dict):
            E, A = properties['E'], properties['A']
        self.type = '2D truss'
        self.tag = tag
        self.E = E
        self.A = A
        self.geo=geo
        self.Q = np.array([0.0])

        self.q = {'1':0}
        self.v = {'1':0}

        self.w = {'1':0.0}

    def __repr__(self):
        return 'truss-{}'.format(self.tag)

    def N(self):
        L = self.L
        N1 = Polynomial([1,-1/L])
        N2 = Polynomial([0, 1/L])
        return np.array([N1,N2])

    def B(self):
        L = self.L
        B1 = self.N()[0].deriv(1)
        B2 = self.N()[1].deriv(1)
        return np.array([B1,B2])

    def v0_vector(self):
        EA  = self.E*self.A
        L   = self.L
        e0  = self.e0
        q0  = self.q0
        w   = self.w
        v0  =  np.array([e0['1']*L])
        v0 += [q0['1']*L/EA]
        v0 += [w['1']*L*L/(2*EA)]
        return v0

    def q0_vector(self):
        EA = self.E*self.A
        e0 = self.e0['1']
        q0 = self.q0['1']
        q0 = q0 - EA*e0
        return [q0]

    def pw_vector(self):
        L = self.L
        pw = self.w['1']
        N = self.N()
        if isinstance(pw,np.polynomial.Polynomial) or isinstance(pw,float):
            p1 = -quad(N[0]*pw,0,L)[0]
            p2 = -quad(N[1]*pw,0,L)[0]
        else:
            print('Unsupported load vector pw')

        return np.array([p1,0,p2,0])

    def ke_matrix(self):
        ag = self.ag()
        k = self.k_matrix()
        return ag.T@(k*ag)

    def kg_matrix(self, N): 
        """return element local stiffness Matrix"""
        E = self.E
        A = self.A
        L = self.L
        k = Structural_Matrix([N/L])
        # Metadata
        k.tag = 'kg'
        k.row_data = ['q_'+ key for key in self.q0.keys()]
        k.column_data = ['v_' + key for key in self.e0.keys()]
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def k_matrix(self): 
        """return element local stiffness Matrix"""

        E = self.E
        A = self.A
        L = self.L
        k = Structural_Matrix([E*A/L])
        # Metadata
        k.tag = 'k'
        k.row_data = ['q_'+ key for key in self.q0.keys()]
        k.column_data = ['v_' + key for key in self.e0.keys()]
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def bg_matrix(self): 
        """return element static matrix, :math:`\\mathbf{b}_g`"""

        cs = self.cs
        sn = self.sn
        bg = np.array([[-cs],
                       [-sn],
                       [ cs],
                       [ sn]])
        return bg

    def f_matrix(self, Roption=True):
        """return element flexibility matrix, :math:`\\mathbf{f}`"""

        A = self.A
        E = self.E
        L = self.L
        f = Structural_Matrix([L/(E*A)])

        if Roption:
            pass

        # Metadata
        f.tag = 'f'
        f.column_data = ['q_' + key for key in self.q0.keys()]
        f.row_data = ['v_'+ key for key in self.e0.keys()]
        f.c_cidx = f.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]

        return f

    def ag(self):
        cs = self.cs
        sn = self.sn
        ag = np.array([[-cs,-sn , cs, sn],])
        return ag


    def GLstrain(self):
        Li = self.Li
        L0 = self.L0
        E_GL = (Li**2 - L0**2) / (2*L0**2)
        return E_GL

    def iGLstrain(self):
        """incremental Green-Lagrange strain"""
        L  = self.L
        Li = self.Li
        E_GL = (L**2 - Li**2) / (2*Li**2)
        return E_GL

class TaperedTruss(Truss):
    def k_matrix(self):
        if isinstance(self.E,float):
            E = Polynomial([self.E])
        else:
            E = self.E
        if isinstance(self.A,float):
            A = Polynomial([self.A])
        else:
            A = self.A
        A = self.A 
        L = self.L
        B = self.B()
        # ke = np.zeros((self.ndf,self.ndf))
        # for i in range(self.ndf):
        #     for j in range(self.ndf):
        #         f = lambda x: B[i](x)*E(x)*A(x)*B[j](x)
        #         ke[i,j] = quad(f,0,L)[0]

        f = lambda x: B[0](x)*E(x)*A(x)*B[0](x)
        k = Structural_Matrix([quad(f,0,L)[0]])

        # Metadata
        k.tag = 'k'
        k.row_data = ['q_'+ key for key in self.q0.keys()]
        k.column_data = ['v_' + key for key in self.e0.keys()]
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def q0_vector(self):
        # EA = self.E*self.A()
        # e0 = self.e0['1']
        # q0 = self.q0['1']
        # q0 = q0 - EA*e0
        return [0]

    def v0_vector(self):
        return [0]

class Truss3D(Element):
    ndf = 3
    ndm = 3
    force_dict = {'1':0}
    def __init__(self, tag, iNode, jNode, A, E):
        super().__init__(self.ndf, self.ndm, self.force_dict, [iNode,jNode])
        self.type = '2D truss'
        self.tag = tag
        self.A = A
        self.E = E

        self.q = {'1':0}
        self.v = {'1':0}

    def __repr__(self):
        return 'tr-{}'.format(self.tag)

    def bg_matrix(self):
        """return element static matrix, bg - pp. 57"""
        cs = self.cs
        sn = self.sn
        cz = self.cz
        bg = np.array([[-cs],
                       [-sn],
                       [-cz],
                       [ cs],
                       [ sn],
                       [ cz]])
        return bg


class Frame(Element):
    """linear 2D Euler-Bernouli frame element"""

    nv = 3
    nn = 2
    nq = 3
    ndm = 2
    ndf = 3
    force_dict = {'1':0, '2':0, '3': 0}
    Qpl = np.zeros((2,nq))

    def __init__(self, tag, iNode, jNode, E=None, A=None, I=None,properties=None):
        super().__init__(self.ndf, self.ndm, self.force_dict, [iNode,jNode])
        if isinstance(properties,dict):
            E, A, I = properties['E'], properties['A'], properties['I']
        self.type = '2D beam'
        self.tag = tag
        self.E = E
        self.A = A
        self.I = I
        self.q = {'1':0, '2':0, '3': 0}     # basic element force
        self.v = {'1':0, '2':0, '3': 0}
        self.w = {'x':0, 'y':0}             #':  "uniform element loads in x and y",

    def __repr__(self):
        return 'el-{}'.format(self.tag) 

    def elastic_curve(self, x, end_rotations, scale=10, global_coord=False):
        n = len(x)
        L = self.L
        N1 = 1-3*(x/L)**2+2*(x/L)**3
        N2 = L*(x/L-2*(x/L)**2+(x/L)**3)
        N3 = 3*(x/L)**2-2*(x/L)**3
        N4 = L*((x/L)**3-(x/L)**2)
        vi = end_rotations[0]
        vj = end_rotations[1]
        y = np.array(vi*N2+vj*N4)*scale
        xy = np.concatenate(([x],[y]))
        if global_coord:
            x0 = self.nodes[0].x
            y0 = self.nodes[0].y
            Rz = self.Rz_matrix()[0:2,0:2]
            xy = Rz@xy + [[x0]*n,[y0]*n]
        return xy

    @property
    def enq(self): 
        """element number of forces, considering deprecation"""
        # output like [1, 1, 1]
        return [1 for x in self.q]

    def ag(self):
        cs = self.cs
        sn = self.sn
        L = self.L
        ag = np.array([[ -cs ,  -sn , 0,  cs ,   sn , 0],
                       [-sn/L,  cs/L, 1, sn/L, -cs/L, 0],
                       [-sn/L,  cs/L, 0, sn/L, -cs/L, 1]])

        if self.dofs[0] == self.dofs[3] or self.dofs[1] == self.dofs[4]:
            ag[0,:] = [0.0]*ag.shape[1]
        return ag

    def ah(self):
        MR = [1 if x else 0 for x in self.rel.values()]
        ah = np.array([[1-MR[0],          0          ,             0        ],
                       [   0   ,       1-MR[1]       ,  -0.5*(1-MR[2])*MR[1]],
                       [   0   , -0.5*(1-MR[1])*MR[2],          1-MR[2]     ]])
        return ah

    def k_matrix(self): 
        """return element local stiffness Matrix"""
        E = self.E
        A = self.A
        I = self.I
        L = self.L
        EI = E*I
        ah = self.ah()
        k = np.array([[E*A/L,    0   ,   0   ],
                      [  0  , 4*EI/L , 2*EI/L],
                      [  0  , 2*EI/L , 4*EI/L]])
        k = ah.T @ k @ ah

        # Assemble matrix metadata
        k = Structural_Matrix(k)
        k.tag = 'k'
        k.row_data = ['q_'+ key for key in self.q0.keys()]
        k.column_data = ['v_' + key for key in self.v0.keys()]
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def f_matrix(self, Roption=False):
        """Flexibility matrix of an element.

        """
        EA = self.E*self.A
        EI = self.E*self.I
        L = self.L
        f = Structural_Matrix([
            [L/EA,     0    ,      0   ],
            [  0 ,  L/(3*EI), -L/(6*EI)],
            [  0 , -L/(6*EI),  L/(3*EI)]])
        ide = set(range(3))

        if Roption:  
            if self.rel["2"]:
                f[:,1] = [0.0]*f.shape[0]
                f[1,:] = [0.0]*f.shape[1]

            if self.rel["3"]:
                f[:,2] = [0.0]*f.shape[0]
                f[2,:] = [0.0]*f.shape[1]

        # Define matrix metadata
        f.tag = 'f'
        f.column_data = ['q_'+ key for key in self.q0.keys()]
        f.row_data = ['v_' + key for key in self.v0.keys()]
        f.c_cidx = f.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return f

    def ke_matrix(self):
        """return element global stiffness Matrix"""

        k  = self.k_matrix()
        ag = self.ag()

        ke = ag.T @ k @ ag
        ke = Structural_Matrix(ke) 
        ke.row_data = ke.column_data = ['u_'+str(int(dof)) for dof in self.dofs]
        return ke

    def bg_matrix(self, Roption=False): 
        """return element global static matrix, bg"""

        cs = self.cs
        sn = self.sn
        L  = self.L
        #                x      ri      rj   Global
        bg = np.array([[-cs, -sn/L,  -sn/L],  # x
                       [-sn,  cs/L,   cs/L],  # y
                       [0.0,   1.0,    0.0],  # rz
                       [ cs,  sn/L,   sn/L],
                       [ sn, -cs/L,  -cs/L],
                       [0.0,   0.0,    1.0]])
        if Roption:
            if self.rel['2']:
                bg[:,1] = [0.0]*bg.shape[0]

            if self.rel['3']:
                bg[:,2] = [0.0]*bg.shape[0]

        if self.dofs[0] == self.dofs[3] or self.dofs[1] == self.dofs[4]:
            bg[:,0] = [0.0]*bg.shape[0]
            # bg = np.delete(bg, 0, axis=1)

        bg = Structural_Matrix(bg)
        bg.tag = 'b'
        bg.column_data = ['x', 'ri', 'rj']
        bg.row_data = ['x', 'y', 'rz']
        bg.c_cidx = bg.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]


        return bg

    def v0_vector(self):
        EA = self.E*self.A
        EI = self.E*self.I
        L = self.L
        e0 = self.e0
        w = self.w
        v0 =  np.array([e0['1']*L, -e0['2']*L/2, e0['2']*L/2])
        v0 += np.array([w['x']*L*L/(2*EA), w['y']*L**3/(24*EI), -w['y']*L**3/(24*EI)])
        v0 = Structural_Matrix(v0)
        v0.tag = 'v'
        v0.column_data = ['v_0']
        v0.row_data = ['axial', 'i-rotation', 'j-rotation']
        v0.c_cidx = False
        v0.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return v0 

    def q0_vector(self):
        L = self.L
        E = self.E
        A = self.A
        I = self.I
        e0= self.e0
        w = self.w
        q0 =  np.array([-e0['1']*E*A, +e0['2']*E*I, -e0['2']*E*I])
        if self.rel['2']:
            q0[1] *= 0
            q0[1] *= 1.5 
        if self.rel['3']:
            q0[2] *= 0
            q0[2] *= 1.5 

        q0[1] += -w['y']*L**2/12*((not self.rel['2']) and (not self.rel['3'])) - w['y']*L**2/8*(self.rel['3'])
        q0[2] += +w['y']*L**2/12*((not self.rel['2']) and (not self.rel['3'])) + w['y']*L**2/8*(self.rel['2'])

        # Metadata
        q0 = Structural_Matrix(q0)
        q0.tag = 'q'
        q0.column_data = ['q_0']
        q0.row_data = ['axial', 'M_i', 'M_j']
        q0.c_cidx = False
        q0.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]

        return q0 

    def κ(self, k, state=None, form=None):
        """Defines element curvature and calculates the resulting end deformations.

        """
        if form ==None: form = 'uniform'
        if state==None: state = -1
        L = self.L
        table = {
            'uniform':   [-0.5*k*L, 0.5*k*L],
            'ilinear':   [-k/3*L, k/6*L],
            'jlinear':   [-k/6*L, k/2*L],
            'parabolic': [-k/3*L, k/3*L],
        }
        self.nodes[0].model.states[state]['v0'][self.tag]['2'] = table[form][0]
        self.nodes[0].model.states[state]['v0'][self.tag]["3"]= table[form][1]
        return


class PrismFrame(Frame):
    axial_force = 0.0 

    def k_matrix_exact(self, axial_force=0.0):
        P = axial_force

        L = self.L
        EI = self.E*self.I
        k = np.sqrt(abs(P/EI))
        kL = k*L
        d = 2*(1-np.cos(kL))-(kL*np.sin(kL))
        ka = self.E*self.A / L
        kfv = (EI/(L**3))*(kL**3)*np.sin(kL)/d
        kmv = (EI/(L**2))*((kL**2)*(1-np.cos(kL))/d)
        kft = kmv
        kmt =(EI/L)*((kL*(np.sin(kL)-kL*np.cos(kL)))/d)
        kmth = (EI/L)*((kL*(kL-np.sin(kL)))/d)

        eltk = np.array([ 
            [  ka,  0.0,   0.0, -ka,   0.0,   0.0],
            [ 0.0,  kfv,   kft,  0.0, -kfv,   kft],
            [ 0.0,  kmv,   kmt,  0.0, -kmv,  kmth],
            [ -ka,  0.0,   0.0,  ka,   0.0,   0.0],
            [ 0.0, -kfv,  -kft,  0.0,  kfv,  -kft],
            [ 0.0,  kft,  kmth,  0.0, -kft,   kmt]])

        # Assemble matrix metadata
        k = Structural_Matrix(eltk)
        k.tag = 'k'
        k.row_data = ['q_'+ key for key in self.q0.keys()]*2
        k.column_data = ['v_' + key for key in self.v0.keys()]*2
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def k_matrix(self, axial_force=None):
        """Evaluates the 2D beam-column stiffness matrix using a
        power series expansion
        """
        if axial_force is None: P = self.axial_force
        else: P = axial_force
        L = self.L
        EI = self.E*self.I
        k = np.sqrt(P/EI)
        kL = k*L
        ka = self.E*self.A / L
        kfv = 12*EI/(L**3)-6*EI/(5*L)*k**2-EI*L/700*k**4
        kmv = 6*EI/L**2 - EI/10*k**2 - EI*L**2/1400*k**4
        kft = kmv
        kmt = 4*EI/L-2*EI*L/15*k**2 - 11*EI/6300*L**3*k**4
        kmth = 2*EI/L+EI*L/30*k**2 + 13*EI*L**3*k**4/12600

        eltk = np.array([
            [  ka,  0.0,  -0.0, -ka,   0.0,   0.0],
            [ 0.0,  kfv,   kft, -0.0, -kfv,   kft],
            [ 0.0,  kmv,   kmt, -0.0, -kmv,  kmth],
            [ -ka,  0.0,  -0.0,  ka,   0.0,   0.0],
            [ 0.0, -kfv,  -kft, -0.0,  kfv,  -kft],
            [ 0.0,  kft,  kmth, -0.0, -kft,   kmt]])

        # hinges
        # cols = rows = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        # f = np.linalg.inv(eltk)

        # for row in rows:
        #     f[row, :] = [0.0]*6
        # for col in cols:
        #     f[:, col] = [0.0]*6
        # k = np.linalg.inv(f)


        # Matrix metadata
        k = Structural_Matrix(eltk)
        k.tag = 'k'
        k.row_data = ['q_'+ key for key in self.q0.keys()]*2
        k.column_data = ['v_' + key for key in self.v0.keys()]*2
        k.c_cidx = k.c_ridx = [int(key)-1 for key in self.rel.keys() if self.rel[key]]
        return k

    def ke_matrix(self, axial_force=None):
        T = self.Rz_matrix()
        k = self.k_matrix(axial_force)
        ke = T.T @ k @ T
        return ke

class Beam3d(Element):
    nv = 6
    nn = 2
    ndm = 3
    ndf = 6
    force_dict = {'1':0, '2':0, '3': 0}

    def __init__(self, tag, iNode, jNode, E, A, Iy, Iz, Gs, Kv):
        super().__init__(self.ndf, self.ndm, self.force_dict, [iNode,jNode])
        self.type = '3D beam'
        self.tag = tag
        self.E  =  E
        self.A  =  A
        self.Gs = Gs
        self.Iy = Iy
        self.Iz = Iz
        self.Kv = Kv
        # self.q = {'1':0, '2':0, '3': 0}     # basic element force
        # self.v = {'1':0, '2':0, '3': 0}

        self.w = {'x':0, 'y':0, 'z': 0}             #':  "uniform element loads in x and y",


    def __repr__(self):
        return 'el-{}'.format(self.tag)

    def k_matrix(self):
        """
        Calculate the stiffness matrix for a 3D elastic Bernoulli
        beam element.

        Parameters
        --------------------------------------------
        E: Young's modulus
        G: Shear modulus
        A: Cross section area
        Iy: Moment of inertia, local y-axis
        Iz: Moment of inertia, local z-axis
        Kv: Saint-Venant's torsion constant

        Returns:
            Kle                      local beam stiffness matrix (12 x 12)
        """

        L,E,Gs,A,Iy,Iz,Kv = self.L, self.E, self.Gs, self.A, self.Iy, self.Iz, self.Kv

        a = E*A/L
        b = 12*E*Iz/L**3
        c = 6*E*Iz/L**2
        d = 12*E*Iy/L**3
        e = 6*E*Iy/L**2
        f = Gs*Kv/L
        g = 2*E*Iy/L
        h = 2*E*Iz/L

        Kle = np.mat([
            [ a,  0,  0,  0,  0,   0,  -a,  0,  0,  0,  0,   0  ],
            [ 0,  b,  0,  0,  0,   c,   0, -b,  0,  0,  0,   c  ],
            [ 0,  0,  d,  0, -e,   0,   0,  0, -d,  0, -e,   0  ],
            [ 0,  0,  0,  f,  0,   0,   0,  0,  0, -f,  0,   0  ],
            [ 0,  0, -e,  0, 2*g,  0,   0,  0,  e,  0,  g,   0  ],
            [ 0,  c,  0,  0,  0,   2*h, 0, -c,  0,  0,  0,   h  ],
            [-a,  0,  0,  0,  0,   0,   a,  0,  0,  0,  0,   0  ],
            [ 0, -b,  0,  0,  0,  -c,   0,  b,  0,  0,  0,  -c  ],
            [ 0,  0, -d,  0,  e,   0,   0,  0,  d,  0,  e,   0  ],
            [ 0,  0,  0, -f,  0,   0,   0,  0,  0,  f,  0,   0  ],
            [ 0,  0, -e,  0,  g,   0,   0,  0,  e,  0, 2*g,  0  ],
            [ 0,  c,  0,  0,  0,   h,   0, -c,  0,  0,  0,  2*h ]])

        # fle = L/2*np.mat([qx, qy, qz, qw, -qz*L/6, qy*L/6, qx, qy, qz, qw, qz*L/6, -qy*L/6]).T

        # n2 = np.array([0.,0.,0.])
        # n2[0] =  n3[1]*n1[2]-n3[2]*n1[1]
        # n2[1] = -n1[2]*n3[0]+n1[0]*n3[2]
        # n2[2] =  n3[0]*n1[1]-n1[0]*n3[1]

        #An = np.append([n1,n2],[n3],0)

        # G = np.mat([
        #     [ n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        #     [ n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        #     [ n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
        #     [ 0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0    ],
        #     [ 0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0    ],
        #     [ 0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0    ],
        #     [ 0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0    ],
        #     [ 0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0    ],
        #     [ 0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0    ],
        #     [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2]],
        #     [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2]],
        #     [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2]]
        # ])

        # Ke = G.T*Kle*G
        # fe = G.T*fle
        return Kle
        # if eq == None:
        #     return Kle
        # else:
        #     return Kle,fe


    def ke_matrix(self):
        Kle = self.k_matrix()
        ex = [self.nodes[0].x, self.nodes[1].x]
        ey = [self.nodes[0].y, self.nodes[1].y]
        ez = [self.nodes[0].z, self.nodes[1].z]
        b = np.mat([
                    [ex[1]-ex[0]],
                    [ey[1]-ey[0]],
                    [ez[1]-ez[0]]])

        L = self.L
        n1 = np.asarray(b.T/L).reshape(3,)

        eo = np.asmatrix(eo)
        lc = np.asscalar(np.sqrt(eo*eo.T))
        n3 = np.asarray(eo/lc).reshape(3,)


        n2 = np.array([0.,0.,0.])
        n2[0] =  n3[1]*n1[2]-n3[2]*n1[1]
        n2[1] = -n1[2]*n3[0]+n1[0]*n3[2]
        n2[2] =  n3[0]*n1[1]-n1[0]*n3[1]

        T = np.mat([
            [ n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
            [ n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
            [ n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0,     0,     0,     0    ],
            [ 0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0,     0,     0,     0    ],
            [ 0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0,     0,     0,     0    ],
            [ 0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0,     0,     0,     0    ],
            [ 0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2], 0,     0,     0    ],
            [ 0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2], 0,     0,     0    ],
            [ 0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2], 0,     0,     0    ],
            [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n1[0], n1[1], n1[2]],
            [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n2[0], n2[1], n2[2]],
            [ 0,     0,     0,     0,     0,     0,     0,     0,     0,     n3[0], n3[1], n3[2]]
        ])
        Ke = T.T*Kle*T
        return Ke


