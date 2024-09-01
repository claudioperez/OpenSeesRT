import opensees.backend as anp

# 1. Expressions
#------------------------------

# constant parameters
dw, tw = 18, 18

def area(tf, bf):
    return bf*tf + dw*tw

def moi(tf, bf):
    area = bf*tf + dw *tw
    y_c = 1/(2*area) * (tw*(dw+tf)**2 + (bf-tw)*tf**2)
    return tw*(dw+tf)**3/3 + (bf - tw)*tf**3/3 - area*y_c**2

# 2. Components/elements
#------------------------------
import opensees as ops
from opensees import element
corot_transf = elle.beam2d.transform_no2(elle.beam2d.geom_no1)
elastic_beam  = elle.beam2d.resp_no1


# Create model Assembler
model = ops.SkeletalModel(ndm=2,ndf=3)

# Add parameters and expressions to the assembler
tf, bf = model.param("tf",  "bf")
A      = model.expr(area, tf, bf)
I      = model.expr(moi,  tf, bf)

# Define basic element templates
basic_girder  = elastic_beam(A=A, E=3600., I=I)
basic_column  = elastic_beam(I=30**4/12, A=30**2, E=3600.0)

girder = corot_transf(basic_girder)
column = corot_transf(basic_column)

# 3. Assembly
#------------------------------
ft = 12
B, H = 30.*ft, 13.*ft
# Set up nodes
model.node("1",  0.,  0.)
model.node("2",  0.,  H )
model.node("3", B/2,  H )
model.node("4",  B ,  H )
model.node("5",  B ,  0.)

#         template      nodes  identifier
model.elem(column, ["1", "2"],  name="a")
model.elem(girder, ["2", "3"],  name="b")
model.elem(girder, ["3", "4"],  name="c")
model.elem(column, ["4", "5"],  name="d")

#         node  fixities
model.boun("1", [1,1,1])
model.boun("5", [1,1,1])

#         node  magnitude  direction
model.load("2",    1000.0, dof="x")
model.load("2",      -2.0, dof="y")
model.load("4",      -2.0, dof="y")

#f = model.compose(_jit_force=True)
#f = ops.analysis.Displacement(_jit_force=True)


