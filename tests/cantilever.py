import opensees
import opensees.library as lib
from opensees.tcl import dumps

from numpy import linspace

L  = 1
ne = 10

transform = lib.LinearTransform(vecxz=[0,0,1])

element = lib.element.ElasticBeamColumn3D(nodes=[1,2],
    E   = 1,
    A   = 1,
    G   = 1,
    J   = 1,
    ixc = 2,
    iyc = 2,
    transform=transform,
)

model = opensees.model(
    ndf=6,
    nodes={i+1: [x, 0, 0] for i,x in enumerate(linspace(0, L, ne+1))},
    elems={i+1: [element, (i+1, i+2)] for i in range(ne)},
    bound={1: [1]*6}
)

# print(dumps(model))

#
# Patterns
#
loads = {
        2: [0.0]*3+[0.0, 1.0, 0.0]
}

patterns = {
        1: opensees.library.pattern.Plain(1, "Linear", loads)
}


#
# Analysis
#
strategy = {
  "constraints": ["Plain"],                   # how it handles boundary conditions
  "numberer":    ["Plain"],                   # renumber dof's to minimize band-width
  "system":      ["BandGeneral"],             # how to store and solve the system of equations in the analysis
  "test":        ["NormDispIncr", 1e-8, 6],   # determine if convergence has been achieved at the end of an iteration step
  "algorithm":   ["Newton"],                  # use Newton's solution algorithm: updates tangent stiffness at every iteration
  "integrator":  ["LoadControl", 0.01],       # determine the next time step for an analysis
}

analysis = opensees.analysis.StaticAnalysis(model, strategy, patterns)

for i in range(10):
    analysis.analyze(1)
    print(analysis.rt.getNodeResponse(2, "displ"))


