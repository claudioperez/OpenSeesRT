import opensees
import opensees.library as g3


transform = g3.LinearTransform(vecxz=[0,0,1])

beam = g3.element.ElasticBeamColumn3D(nodes=[1,2],
    E   = 1,
    A   = 1,
    G   = 1,
    J   = 1,
    ixc = 2,
    iyc = 2,
    transform=transform,
)

nodes = {
     1:  [0.0, 0, 0],
     2:  [0.1, 0, 0],
     3:  [0.2, 0, 0],
     4:  [0.3, 0, 0],
     5:  [0.4, 0, 0],
     6:  [0.5, 0, 0],
     7:  [0.6, 0, 0],
     8:  [0.7, 0, 0],
     9:  [0.8, 0, 0],
     10: [0.9, 0, 0],
     11: [1.0, 0, 0]
}

elems = {
     1:  [beam, ( 1,  2)],
     2:  [beam, ( 2,  3)],
     3:  [beam, ( 3,  4)],
     4:  [beam, ( 4,  5)],
     5:  [beam, ( 5,  6)],
     6:  [beam, ( 6,  7)],
     7:  [beam, ( 7,  8)],
     8:  [beam, ( 8,  9)],
     9:  [beam, ( 9, 10)],
     10: [beam, (10, 11)]
}

model = opensees.model(ndf=6, nodes=nodes, elems=elems, bound={1: [1]*6})

tcl_script = opensees.tcl.dumps(model)


#-------------------------------------------------------------------------
# Tcl Runtime
#-------------------------------------------------------------------------

print(tcl_script)

tcl = opensees.tcl.TclRuntime()

tcl.eval(tcl_script)

tcl.eval("print -json")



#-------------------------------------------------------------------------
# Analysis
#-------------------------------------------------------------------------
#
# Loads & Patterns
#
loads    = {
    2: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
}

patterns = {
    # TODO: remove tag from Plain?
    1: g3.pattern.Plain(1, "Linear", loads)
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

# analysis = opensees.analysis.StaticAnalysis(model, strategy, patterns)

# analysis.analyze(3)

# print(analysis.rt.getNodeResponse(2, "displ"))

