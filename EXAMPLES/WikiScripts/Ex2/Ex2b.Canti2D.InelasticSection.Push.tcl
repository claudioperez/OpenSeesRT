# --------------------------------------------------------------------------------------------------
# Example 2. 2D cantilever column, static pushover
# element
#			Silvia Mazzoni & Frank McKenna, 2006
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
wipe;					# clear memory of all past model definitions
file mkdir Data; 				# create data directory
model BasicBuilder -ndm 2 -ndf 3;		# Define the model builder, ndm=#dimension, ndf=#dofs


# define GEOMETRY -------------------------------------------------------------
set LCol 432; 		# column length
set Weight 2000.; 		# superstructure weight
# define section geometry
set HCol 60; 		# Column Depth
set BCol 60;		# Column Width

# calculated parameters
set PCol $Weight; 		# nodal dead-load weight per column
set g 386.4;			# g.
set Mass [expr $PCol/$g];		# nodal mass
# calculated geometry parameters
set ACol [expr $BCol*$HCol*1000];				# cross-sectional area, make stiff
set IzCol [expr 1./12.*$BCol*pow($HCol,3)]; 			# Column moment of inertia

# nodal coordinates:
node 1 0 0;			# node#, X, Y
node 2 0 $LCol 		

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			# node DX DY RZ

# nodal masses:
mass 2 $Mass  1e-9 0.;		# node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

# Define ELEMENTS & SECTIONS  -------------------------------------------------------------
set ColMatTagFlex 2;			# assign a tag number to the column flexural behavior
set ColMatTagAxial 3;			# assign a tag number to the column axial behavior	
set ColSecTag 1;				# assign a tag number to the column section tag
set BeamSecTag 2;				# assign a tag number to the beam section tag

# MATERIAL parameters
set fc -4.; 				# CONCRETE Compressive Strength (+Tension, -Compression)
set Ec [expr 57*sqrt(-$fc*1000)]; 		# Concrete Elastic Modulus (the term in sqr root needs to be in psi

# COLUMN section
# calculated stiffness parameters
set EICol [expr $Ec*$IzCol];			# EI, for moment-curvature relationship
set EACol [expr $Ec*$ACol];			# EA, for axial-force-strain relationship
set MyCol 130000;				# yield moment
set PhiYCol 0.65e-4;			# yield curvature
set EIColCrack [expr $MyCol/$PhiYCol];	# cracked section inertia
set b 0.01 ;				# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $ColMatTagFlex $MyCol $EIColCrack $b; 		# bilinear behavior for flexure
uniaxialMaterial Elastic $ColMatTagAxial $EACol;				# this is not used as a material, this is an axial-force-strain response
section Aggregator $ColSecTag $ColMatTagAxial P $ColMatTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set ColTransfTag 1; 			# associate a tag to column transformation
geomTransf Linear $ColTransfTag  ; 	

# element connectivity:
set numIntgrPts 5;								# number of integration points for force-based element
element nonlinearBeamColumn 1 1 2 $numIntgrPts $ColSecTag $ColTransfTag;	# self-explanatory when using variables

# Define RECORDERS -------------------------------------------------------------
recorder Node -file Data/DFree.out -time -node 2 -dof 1 2 3 disp;		# displacements of free nodes
recorder Node -file Data/DBase.out -time -node 1 -dof 1 2 3 disp;		# displacements of support nodes
recorder Node -file Data/RBase.out -time -node 1 -dof 1 2 3 reaction;		# support reaction
recorder Drift -file Data/Drift.out -time -iNode 1 -jNode 2 -dof 1   -perpDirn 2 ;	# lateral drift
recorder Element -file Data/FCol.out -time -ele 2 globalForce;						# element forces -- column
recorder Element -file Data/ForceColSec1.out -time -ele 1 section 1 force;				# Column section forces, axial and moment, node i
recorder Element -file Data/DefoColSec1.out -time -ele 1 section 1 deformation;				# section deformations, axial and curvature, node i
recorder Element -file Data/ForceColSec$numIntgrPts.out -time -ele 1 section $numIntgrPts force;		# section forces, axial and moment, node j
recorder Element -file Data/DefoColSec$numIntgrPts.out -time -ele 1 section $numIntgrPts deformation;		# section deformations, axial and curvature, node j


# define GRAVITY -------------------------------------------------------------
pattern Plain 1 Linear {
   load 2 0 -$PCol 0
}

# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;			# convergence tolerance for test
constraints Plain;     		# how it handles boundary conditions
numberer Plain;			# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;		# how to store and solve the system of equations in the analysis
test NormDispIncr $Tol 6 ; 		# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;			# use Newton's solution algorithm: updates tangent stiffness at every iteration
set NstepGravity 10;  		# apply gravity in 10 steps
set DGravity [expr 1./$NstepGravity]; 	# first load increment;
integrator LoadControl $DGravity;	# determine the next time step for an analysis
analysis Static;			# define type of analysis static or transient
analyze $NstepGravity;		# apply gravity

# ------------------------------------------------- maintain constant gravity loads and reset time to zero
loadConst -time 0.0

puts "Model Built"

# STATIC PUSHOVER ANALYSIS --------------------------------------------------------------------------------------------------
#
# we need to set up parameters that are particular to the model.
set IDctrlNode 2;			# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement contro
set Dmax [expr 0.05*$LCol];		# maximum displacement of pushover. push to 10% drift.
set Dincr [expr 0.001*$LCol];		# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis

# create load pattern for lateral pushover load
set Hload $Weight;				# define the lateral load as a proportion of the weight so that the pseudo time equals the lateral-load coefficient when using linear load pattern
pattern Plain 200 Linear {;			# define load pattern -- generalized
	load 2 $Hload 0.0 0.0 0.0 0.0 0.0;	# define lateral load in static lateral analysis
}

# ----------- set up analysis parameters
# CONSTRAINTS handler -- Determines how the constraint equations are enforced in the analysis (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/617.htm)
#          Plain Constraints -- Removes constrained degrees of freedom from the system of equations (only for homogeneous equations)
#          Lagrange Multipliers -- Uses the method of Lagrange multipliers to enforce constraints 
#          Penalty Method -- Uses penalty numbers to enforce constraints --good for static analysis with non-homogeneous eqns (rigidDiaphragm)
#          Transformation Method -- Performs a condensation of constrained degrees of freedom 
constraints Plain;		

# DOF NUMBERER (number the degrees of freedom in the domain): (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/366.htm)
#   determines the mapping between equation numbers and degrees-of-freedom
#          Plain -- Uses the numbering provided by the user 
#          RCM -- Renumbers the DOF to minimize the matrix band-width using the Reverse Cuthill-McKee algorithm 
numberer Plain

# SYSTEM (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/371.htm)
#   Linear Equation Solvers (how to store and solve the system of equations in the analysis)
#   -- provide the solution of the linear system of equations Ku = P. Each solver is tailored to a specific matrix topology. 
#          ProfileSPD -- Direct profile solver for symmetric positive definite matrices 
#          BandGeneral -- Direct solver for banded unsymmetric matrices 
#          BandSPD -- Direct solver for banded symmetric positive definite matrices 
#          SparseGeneral -- Direct solver for unsymmetric sparse matrices 
#          SparseSPD -- Direct solver for symmetric sparse matrices 
#          UmfPack -- Direct UmfPack solver for unsymmetric matrices 
system BandGeneral

# TEST: # convergence test to 
# Convergence TEST (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/360.htm)
#   -- Accept the current state of the domain as being on the converged solution path 
#   -- determine if convergence has been achieved at the end of an iteration step
#          NormUnbalance -- Specifies a tolerance on the norm of the unbalanced load at the current iteration 
#          NormDispIncr -- Specifies a tolerance on the norm of the displacement increments at the current iteration 
#          EnergyIncr-- Specifies a tolerance on the inner product of the unbalanced load and displacement increments at the current iteration 
set Tol 1.e-8;                        # Convergence Test: tolerance
set maxNumIter 6;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlag 0;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set TestType EnergyIncr ;	# Convergence-test type
test $TestType $Tol $maxNumIter $printFlag;

# Solution ALGORITHM: -- Iterate from the last time step to the current (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/682.htm)
#          Linear -- Uses the solution at the first iteration and continues 
#          Newton -- Uses the tangent at the current iteration to iterate to convergence 
#          ModifiedNewton -- Uses the tangent at the first iteration to iterate to convergence 
set algorithmType Newton
algorithm $algorithmType;        

# Static INTEGRATOR: -- determine the next time step for an analysis  (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/689.htm)
#          LoadControl -- Specifies the incremental load factor to be applied to the loads in the domain 
#          DisplacementControl -- Specifies the incremental displacement at a specified DOF in the domain 
#          Minimum Unbalanced Displacement Norm -- Specifies the incremental load factor such that the residual displacement norm in minimized 
#          Arc Length -- Specifies the incremental arc-length of the load-displacement path 
# Transient INTEGRATOR: -- determine the next time step for an analysis including inertial effects 
#          Newmark -- The two parameter time-stepping method developed by Newmark 
#          HHT -- The three parameter Hilbert-Hughes-Taylor time-stepping method 
#          Central Difference -- Approximates velocity and acceleration by centered finite differences of displacement 
integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr

# ANALYSIS  -- defines what type of analysis is to be performed (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/324.htm)
#          Static Analysis -- solves the KU=R problem, without the mass or damping matrices. 
#          Transient Analysis -- solves the time-dependent analysis. The time step in this type of analysis is constant. The time step in the output is also constant. 
#          variableTransient Analysis -- performs the same analysis type as the Transient Analysis object. The time step, however, is variable. This method is used when 
#                 there are convergence problems with the Transient Analysis object at a peak or when the time step is too small. The time step in the output is also variable.
analysis Static

#  ---------------------------------    perform Static Pushover Analysis
set Nsteps [expr int($Dmax/$Dincr)];        # number of pushover analysis steps
set ok [analyze $Nsteps];                # this will return zero if no convergence problems were encountered

# ---------------------------------- in case of convergence problems
if {$ok != 0} {      
# change some analysis parameters to achieve convergence
# performance is slower inside this loop
	set ok 0;
	set controlDisp 0.0;		# start from zero
	set D0 0.0;		# start from zero
	set Dstep [expr ($controlDisp-$D0)/($Dmax-$D0)]
	while {$Dstep < 1.0 && $ok == 0} {	
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr ($controlDisp-$D0)/($Dmax-$D0)]
		set ok [analyze 1 ]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 2000  0
			algorithm Newton -initial
			set ok [analyze 1 ]
			test $TestType $Tol $maxNumIter  0
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 ]
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			set ok [analyze 1 ]
			algorithm $algorithmType
		}
	}
	};      # end if ok !0

puts "DonePushover"

