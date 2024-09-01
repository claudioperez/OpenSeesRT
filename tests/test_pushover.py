
script = """
    # -----------------------------------------------------------------------------
    # Example 2. 2D cantilever column, static pushover
    # element

    model basic -ndm 2 -ndf 3;   # Define the model builder, ndm=#dimension, ndf=#dofs


    # define GEOMETRY -------------------------------------------------------------
    set LCol 432;                 # column length
    set Weight 2000.;             # superstructure weight
    # define section geometry
    set HCol 60;                  # Column Depth
    set BCol 60;                  # Column Width

    # calculated parameters
    set PCol $Weight;             # nodal dead-load weight per column
    set g 386.4;                  # g.
    set Mass [expr $PCol/$g];     # nodal mass

    # calculated geometry parameters
    set ACol [expr $BCol*$HCol*1000];                                # cross-sectional area, make stiff
    set IzCol [expr 1./12.*$BCol*pow($HCol,3)];                         # Column moment of inertia

    # nodal coordinates:
    node 1 0 0;                        # node#, X, Y
    node 2 0 $LCol

    # Single point constraints -- Boundary Conditions
    fix 1 1 1 1;                         # node DX DY RZ

    # nodal masses:
    mass 2 $Mass  1e-9 0.;               # node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

    # Define ELEMENTS & SECTIONS  -------------------------------------------------------------
    set ColMatTagFlex 2;                        # assign a tag number to the column flexural behavior
    set ColMatTagAxial 3;                       # assign a tag number to the column axial behavior
    set ColSecTag 1;                            # assign a tag number to the column section tag
    set BeamSecTag 2;                           # assign a tag number to the beam section tag

    # MATERIAL parameters
    set fc -4.;                                 # CONCRETE Compressive Strength (+Tension, -Compression)
    set Ec [expr 57*sqrt(-$fc*1000)];           # Concrete Elastic Modulus (the term in sqr root needs to be in psi

    # COLUMN section
    # calculated stiffness parameters
    set EICol [expr $Ec*$IzCol];                # EI, for moment-curvature relationship
    set EACol [expr $Ec*$ACol];                 # EA, for axial-force-strain relationship
    set MyCol 130000;                           # yield moment
    set PhiYCol 0.65e-4;                        # yield curvature
    set EIColCrack [expr $MyCol/$PhiYCol];      # cracked section inertia
    set b 0.01 ;                                # strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)

    uniaxialMaterial Steel01 $ColMatTagFlex $MyCol $EIColCrack $b;                 # bilinear behavior for flexure
    uniaxialMaterial Elastic $ColMatTagAxial $EACol;                                # this is not used as a material, this is an axial-force-strain response
    section Aggregator $ColSecTag $ColMatTagAxial P $ColMatTagFlex Mz;        # combine axial and flexural behavior into one section (no P-M interaction here)

    # define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
    set ColTransfTag 1;                         # associate a tag to column transformation
    geomTransf Linear $ColTransfTag  ;

    # element connectivity:
    set numIntgrPts 5;                                                                # number of integration points for force-based element
    element nonlinearBeamColumn 1 1 2 $numIntgrPts $ColSecTag $ColTransfTag;        # self-explanatory when using variables

    # Define RECORDERS -------------------------------------------------------------
    recorder Node -file Data/DFree.out -time -node 2 -dof 1 2 3 disp;                # displacements of free nodes
    recorder Node -file Data/DBase.out -time -node 1 -dof 1 2 3 disp;                # displacements of support nodes
    recorder Node -file Data/RBase.out -time -node 1 -dof 1 2 3 reaction;                # support reaction
    recorder Drift -file Data/Drift.out -time -iNode 1 -jNode 2 -dof 1   -perpDirn 2 ;        # lateral drift
    recorder Element -file Data/FCol.out -time -ele 2 globalForce;                                                # element forces -- column
    recorder Element -file Data/ForceColSec1.out -time -ele 1 section 1 force;                                # Column section forces, axial and moment, node i
    recorder Element -file Data/DefoColSec1.out -time -ele 1 section 1 deformation;                                # section deformations, axial and curvature, node i
    recorder Element -file Data/ForceColSec$numIntgrPts.out -time -ele 1 section $numIntgrPts force;                # section forces, axial and moment, node j
    recorder Element -file Data/DefoColSec$numIntgrPts.out -time -ele 1 section $numIntgrPts deformation;                # section deformations, axial and curvature, node j


    # define GRAVITY -------------------------------------------------------------
    pattern Plain 1 Linear {
       load 2 0 -$PCol 0
    }
"""

from opensees import tcl, analysis

rt = tcl.TclRuntime()
rt.eval(script)

NstepGravity = 10;                    # apply gravity, in, 10, steps
DGravity = 1./NstepGravity;           # first load, increment;
Tol      = 1.0e-8;                         # convergence, tolerance, for, test

grav_strategy = {
  "constraints": ["Plain"],                  # how, it, handles, boundary, conditions
  "numberer":    ["Plain"],                  # renumber, dofs, to, minimize, band-width (optimization), if, you, want, to
  "system":      ["BandGeneral"],            # how to store and, solve, the, system, of, equations, in, the, analysis
  "test":        ["NormDispIncr", Tol, 6 ],  # determine, if, convergence, has, been, achieved, at, the, end, of, an, iteration, step
  "algorithm":   ["Newton"],                 # use Newtons, solution, algorithm: updates, tangent, stiffness, at, every, iteration
  "integrator":  ["LoadControl", DGravity],  # determine, the, next, time, step, for, an, analysis
}

grav_analysis = analysis.StaticAnalysis(rt, grav_strategy)

#from opensees.OpenSeesPyRT import libOpenSeesRT
#print(libOpenSeesRT.get_domain(libOpenSeesRT.getRuntime(grav_analysis._rt._interp.tk.interpaddr())))
grav_analysis.analyze(1)                # apply, gravity
print(rt.getNodeResponse(2, "displ"))


pushover_script = """
    #
    # we need to set up parameters that are particular to the model.
    set IDctrlNode 2;			# node where displacement is read for displacement control
    set IDctrlDOF 1;			# degree of freedom of displacement read for displacement contro
    set Dmax [expr 0.05*$LCol];		# maximum displacement of pushover. push to 10% drift.
    set Dincr [expr 0.001*$LCol];		# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis

    # create load pattern for lateral pushover load
    set Hload $Weight;				# define the lateral load as a proportion of the weight so that the pseudo time equals the lateral-load coefficient when using linear load pattern
    pattern Plain 200 Linear { 			# define load pattern -- generalized
    	load 2 $Hload 0.0 0.0 ;         # define lateral load in static lateral analysis
    }
"""

rt.eval(pushover_script)

strategy = {
    "constraints": ["Plain"],
    "numberer"   : ["Plain"],
    "system"     : ["BandGeneral"],
    "test"       : ["EnergyIncr", 1e-8, 6, 0],
    "algorithm"  : ["Newton"],
    "integrator" : ["DisplacementControl",  2, 1, 0.001*10],
}

push_analysis = analysis.StaticAnalysis(rt, strategy)

push_analysis.analyze(3)

print(rt.getNodeResponse(2, "displ"))

