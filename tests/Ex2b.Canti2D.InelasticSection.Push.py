from math import cos,sin,sqrt,pi
import opensees as ops
# --------------------------------------------------------------------------------------------------
# Example 2. 2D cantilever column, static pushover
# element
#                     Silvia Mazzoni & Frank McKenna, 2006
#
#    ^Y
#    |
#    2       __
#    |          |
#    |          |
#    |          |
#  (1)       LCol
#    |          |
#    |          |
#    |          |
#  =1=      _|_  -------->X
#

# SET UP ----------------------------------------------------------------------------
# units: kip, inch, sec
opensees.model('BasicBuilder', ndm=2, ndf=3);              # Define, 'the', model, builder, ndm=#dimension, ndf=#dofs


# define GEOMETRY -------------------------------------------------------------
LCol = 432;               # column, length
Weight = 2000.;               # superstructure, weight
# define section geometry
HCol = 60;               # Column, Depth
BCol = 60;              # Column, Width

# calculated parameters
PCol = Weight;              # nodal, dead-load, weight, per, column
g    = 386.4;               # g.
Mass = PCol/g;              # nodal, mass
# calculated geometry parameters
ACol  = BCol*HCol*1000;                            # cross-sectional, area, make, stiff
IzCol = 1./12.*BCol*HCol**3;                      # Column, moment, of, inertia

# nodal coordinates:
model.node(1, 0 0);                     # node#, X, Y
model.node(2, 0 LCol               )

# Single point constraints -- Boundary Conditions
model.fix(1, 1 1, 1)                      # model.node('DX', DY, RZ)

# nodal masses:
model.mass(2, Mass, 1e-9, 0.);              # node#, Mx, 'My', Mz, Mass=Weight/g, neglect, 'rotational', inertia, 'at', nodes

# Define ELEMENTS & SECTIONS  -------------------------------------------------------------
# ColMatTagFlex = 2;                     # assign, a tag, number, to, the, column, flexural, behavior
# ColMatTagAxial = 3;                     # assign, a tag, number, to, the, column, axial, behavior
# ColSecTag = 1;                            # assign, a tag, number, to, the, column, section, tag

# MATERIAL parameters
fc = -4.;                             # CONCRETE, Compressive, Strength (+Tension, -Compression)
Ec = 57*sqrt(-fc*1000);               # Concrete, Elastic, Modulus (the, term, in, sqr, root, needs, to, be, in, psi

# COLUMN section
# calculated stiffness parameters
EICol = Ec*IzCol;               # EI for moment-curvature relationship
EACol = Ec*ACol;                # EA for axial-force-strain relationship
MyCol = 130000;                 # yield moment
PhiYCol = 0.65e-4;              # yield curvature
EIColCrack = MyCol/PhiYCol;     # cracked section inertia
b = 0.01 ;                      # strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
ColMatFlex  = uniaxial.Steel01(None, MyCol, EIColCrack, b);               # bilinear, 'behavior', for, flexure
ColMatAxial = uniaxial.Elastic(None, EACol);                            # this, 'is', not, 'used', as, a material, this, 'is', an, axial-force-strain, response
ColSec = section.Aggregator(ColSecTag, ColMatAxial, 'P', ColMatFlex, 'Mz');       # combine, 'axial', and, 'flexural', behavior, 'into', one, 'section', (no, P-M, 'interaction', here)

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
ColTransf = opensees.geomTransf.Linear();

# element connectivity:
numIntgrPts = 5;                                                        # number, of, integration, points, for, force-based, element
opensees.element.ForceBeamColumn(1, [1, 2], numIntgrPts, ColSec, ColTransf);       # self-explanatory, 'when', using, variables

# Define RECORDERS -------------------------------------------------------------
# recorder, 'Node', -file, 'Data/DFree.out', -time -model.node(2 -dof, 1 2, 3 disp);              # displacements, 'of', free, nodes
# recorder, 'Node', -file, 'Data/DBase.out', -time -model.node(1 -dof, 1 2, 3 disp);              # displacements, 'of', support, nodes
# recorder, 'Node', -file, 'Data/RBase.out', -time -model.node(1 -dof, 1 2, 3 reaction);              # support, reaction
# recorder, 'Drift', -file, 'Data/Drift.out', -time -iNode, 1 -jNode, 2 -dof, 1   -perpDirn, 2 ;       # lateral, drift
# recorder, 'Element', -file, 'Data/FCol.out', -time -ele, 2 globalForce;                                          # model.element('forces', -- column)
# recorder, 'Element', -file, 'Data/ForceColSec1.out', -time -ele, 1 section, 1 force;                            # Column, 'section', forces, axial, 'and', moment, model.node(i)
# recorder, 'Element', -file, 'Data/DefoColSec1.out', -time -ele, 1 section, 1 deformation;                            # section, deformations, axial, 'and', curvature, model.node(i)
# recorder, 'Element', -file, Data/ForceColSecnumIntgrPts.out -time -ele, 1 section, numIntgrPts, force;              # section, forces, axial, 'and', moment, model.node(j)
# recorder, 'Element', -file, Data/DefoColSecnumIntgrPts.out -time -ele, 1 section, numIntgrPts, deformation;              # section, deformations, axial, 'and', curvature, model.node(j)


# define GRAVITY -------------------------------------------------------------
pattern.Plain(1, 'Linear', {
   load 2 0 PCol=0,
})

# Gravity-analysis parameters -- load-controlled static analysis
NstepGravity = 10;                    # apply gravity, in, 10, steps
DGravity = 1./NstepGravity;           # first load, increment;
Tol      = 1.0e-8;                         # convergence, tolerance, for, test

grav_strategy = {
  "constraints": ["Plain"],                  # how, 'it', handles, 'boundary', conditions
  "numberer":    ["Plain"],                  # renumber, dof's, 'to', minimize, band-width (optimization), if, 'you', want, to
  "system":      ["BandGeneral"],            # how, 'to', store, 'and', solve, 'the', system, 'of', equations, 'in', the, analysis
  "test":        ["NormDispIncr", Tol, 6 ],  # determine, 'if', convergence, 'has', been, 'achieved', at, 'the', end, 'of', an, 'iteration', step
  "algorithm":   ["Newton"],                 # use, Newton's, 'solution', algorithm: updates, 'tangent', stiffness, 'at', every, iteration
  "integrator":  ['LoadControl', DGravity],  # determine, 'the', next, 'time', step, 'for', an, analysis
}

grav_analysis = opensees.analysis.StaticAnalysis(model, grav_strategy)
grav_analysis.analyze(NstepGravity)                # apply, gravity

# ------------------------------------------------- maintain constant gravity loads and retime = to zero
analysis.loadConst(time=0.0)

print("Model Built")

# STATIC PUSHOVER ANALYSIS --------------------------------------------------------------------------------------------------
#
# we need to up = parameters that are particular to the model.
IDctrlNode = 2;                     # model.node(where, displacement, is, read, for, displacement, control)
IDctrlDOF = 1;                     # degree, of, freedom, of, displacement, read, for, displacement, contro
Dmax = 0.05*LCol;              # maximum, displacement, of, pushover., push, to, 10% drift.
Dincr = 0.001*LCol;              # displacement, increment, for, pushover., you, want, this, to, be, very, small, but, not, too, small, to, slow, down, the, analysis

# create load pattern for lateral pushover load
Hload = Weight;                            # define, the, lateral, load, as, a proportion, of, the, weight, so, that, the, pseudo, time, equals, the, lateral-load, coefficient, when, using, linear, load, pattern
pattern.Plain(200, Linear, {                      # define, 'load', pattern -- generalized
    (load 2, Hload, 0.0, 0.0, 0.0, 0.0, 0.0),       # define 'lateral', load 'in', static 'lateral', analysis
})

strategy = {
    "constraints": ['Plain'],
    "numberer"   : ["Plain"],
    "system"     : ["BandGeneral"],
    "test"       : ["EnergyIncr", 1e-8, 6, 0],
    "algorithm"  : ["Newton"],
    "integrator" : ['DisplacementControl',  IDctrlNode, IDctrlDOF, Dincr],
}

analysis = opensees.analysis.StaticAnalysis(model, strategy)

#  ---------------------------------    perform Static Pushover Analysis
Nsteps = int(Dmax/Dincr);        # number, of, pushover, analysis, steps
ok =  analysis.analyze(Nsteps);  # this, will, return, zero, if, no, convergence, problems, were, encountered

# ---------------------------------- in case of convergence problems
if  ok != 0  :
# change some analysis parameters to achieve convergence
# performance is slower inside this loop
       ok = 0;
       controlDisp = 0.0;              # start from zero
       D0 = 0.0;              # start from zero
       Dstep = ((controlDisp-D0)/(Dmax-D0))
       while Dstep < 1.0  and  ok == 0 :
           #controlDisp = [nodeDisp IDctrlNode IDctrlDOF ]
              Dstep =, ((controlDisp-D0)/(Dmax-D0))
              ok = [analysis.analyze(1 ])
              if  ok != 0  :
                     print("Trying Newton with Initial Tangent ..")
                     analysis.test('NormDispIncr',   Tol, 2000, 0)
                     analysis.algorithm('Newton', '-initial')
                     ok = analysis.analyze(1)
                     analysis.test(TestType Tol maxNumIter  0)
                     analysis.algorithm(algorithmType)

              if  ok != 0  :
                     print("Trying Broyden ..")
                     analysis.algorithm('Broyden', 8)
                     ok = analysis.analyze(1)
                     analysis.algorithm(algorithmType)

              if ok != 0  :
                 print("Trying NewtonWithLineSearch ..")
                 analysis.algorithm('NewtonLineSearch', .8)
                 ok = analysis.analyze(1)
                 analysis.algorithm(algorithmType)

print("DonePushover")

