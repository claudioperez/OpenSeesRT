# --------------------------------------------------------------------------------------------------
# Example4. 2D Portal Frame--  Build Model
# nonlinearBeamColumn element, inelastic fiber section
#		Silvia Mazzoni & Frank McKenna, 2006
#
#    ^Y
#    |
#    3_________(3)________4       __ 
#    |                                    |          | 
#    |                                    |          |
#    |                                    |          |
#  (1)                                 (2)       LCol
#    |                                    |          |
#    |                                    |          |
#    |                                    |          |
#  =1=                               =2=      _|_  -------->X
#    |----------LBeam------------|
#

# SET UP ----------------------------------------------------------------------------
wipe;					# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;		# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;				# set up name of data directory
file mkdir $dataDir; 				# create data directory
set GMdir "GMfiles";				# ground-motion file directory
source LibUnits.tcl;				# define basic and system units

# define GEOMETRY -------------------------------------------------------------
set LCol [expr 36*$ft]; 		# column length
set LBeam [expr 42*$ft];		# beam length
set Weight [expr 2000.*$kip]; 		# superstructure weight
# define section geometry
set HCol [expr 5.*$ft]; 		# Column Depth
set BCol [expr 5.*$ft];		# Column Width
set HBeam [expr 8.*$ft];		# Beam Depth
set BBeam [expr 5.*$ft];		# Beam Width

# calculated parameters
set PCol [expr $Weight/2]; 		# nodal dead-load weight per column
set Mass [expr $PCol/$g];		# nodal mass
set MCol [expr 1./12.*($Weight/$LBeam)*pow($LBeam,2)];	# beam-end moment due to distributed load.
# calculated geometry parameters
set ACol [expr $BCol*$HCol];					# cross-sectional area
set ABeam [expr $BBeam*$HBeam];
set IzCol [expr 1./12.*$BCol*pow($HCol,3)]; 			# Column moment of inertia
set IzBeam [expr 1./12.*$BBeam*pow($HBeam,3)]; 		# Beam moment of inertia

# nodal coordinates:
node 1 0 0;			# node#, X, Y
node 2 $LBeam 0
node 3 0 $LCol 		
node 4 $LBeam $LCol 	

# Single point constraints -- Boundary Conditions
fix 1 1 1 0; 			# node DX DY RZ
fix 2 1 1 0; 			# node DX DY RZ
fix 3 0 0 0
fix 4 0 0 0

# nodal masses:
mass 3 $Mass  0. 0.;		# node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes
mass 4 $Mass  0.  0.

# Define ELEMENTS & SECTIONS -------------------------------------------------------------
set ColSecTag 1;			# assign a tag number to the column section	
set BeamSecTag 2;			# assign a tag number to the beam section	
# define section geometry
set coverCol [expr 6.*$in];		# Column cover to reinforcing steel NA.
set numBarsCol 10;			# number of longitudinal-reinforcement bars in each side of column section. (symmetric top & bot)
set barAreaCol [expr 2.25*$in2];	# area of longitudinal-reinforcement bars

# MATERIAL parameters -------------------------------------------------------------------
set IDconcU 1; 			# material ID tag -- unconfined cover concrete
set IDreinf 2; 				# material ID tag -- reinforcement
# nominal concrete compressive strength
set fc [expr -4.0*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec [expr 57*$ksi*sqrt(-$fc/$psi)];	# Concrete Elastic Modulus
# unconfined concrete
set fc1U $fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U -0.003;			# strain at maximum strength of unconfined concrete
set fc2U [expr 0.2*$fc1U];		# ultimate stress
set eps2U	-0.05;			# strain at ultimate stress
set lambda 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
set ftU [expr -0.14*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness
# -----------
set Fy [expr 66.8*$ksi];		# STEEL yield stress
set Es [expr 29000.*$ksi];		# modulus of steel
set Bs 0.01;			# strain-hardening ratio 
set R0 18;				# control the transition from elastic to plastic branches
set cR1 0.925;				# control the transition from elastic to plastic branches
set cR2 0.15;				# control the transition from elastic to plastic branches
uniaxialMaterial Concrete02 $IDconcU $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# build cover concrete (unconfined)
uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;			# build reinforcement material

# FIBER SECTION properties -------------------------------------------------------------
# symmetric section
#                        y
#                        ^
#                        |     
#             ---------------------     --   --
#             |   o     o     o    |     |    -- cover
#             |                       |     |
#             |                       |     |
#    z <--- |          +           |     H
#             |                       |     |
#             |                       |     |
#             |   o     o     o    |     |    -- cover
#             ---------------------     --   --
#             |-------- B --------|
#
# RC section: 
   set coverY [expr $HCol/2.0];	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coverZ [expr $BCol/2.0];	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coreY [expr $coverY-$coverCol]
   set coreZ [expr $coverZ-$coverCol]
   set nfY 16;			# number of fibers for concrete in y-direction
   set nfZ 4;			# number of fibers for concrete in z-direction
   section fiberSec $ColSecTag   {;	# Define the fiber section
	patch quadr $IDconcU $nfZ $nfY -$coverY $coverZ -$coverY -$coverZ $coverY -$coverZ $coverY $coverZ; 	# Define the concrete patch
	layer straight $IDreinf $numBarsCol $barAreaCol -$coreY $coreZ -$coreY -$coreZ;	# top layer reinfocement
	layer straight $IDreinf $numBarsCol $barAreaCol  $coreY $coreZ  $coreY -$coreZ;	# bottom layer reinforcement
    };	# end of fibersection definition

# BEAM section:
section Elastic $BeamSecTag   $Ec $ABeam $IzBeam;	# elastic beam section

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set ColTransfTag 1; 			# associate a tag to column transformation
set BeamTransfTag 2; 			# associate a tag to beam transformation (good practice to keep col and beam separate)
set ColTransfType Linear ;			# options, Linear PDelta Corotational 
geomTransf $ColTransfType $ColTransfTag ; 	# only columns can have PDelta effects (gravity effects)
geomTransf Linear $BeamTransfTag  ; 	

# element connectivity:
set numIntgrPts 5;								# number of integration points for force-based element
element nonlinearBeamColumn 1 1 3 $numIntgrPts $ColSecTag $ColTransfTag;	# self-explanatory when using variables
element nonlinearBeamColumn 2 2 4 $numIntgrPts $ColSecTag $ColTransfTag;
element nonlinearBeamColumn 3 3 4 $numIntgrPts $BeamSecTag $BeamTransfTag;

# Define RECORDERS -------------------------------------------------------------
recorder Node -file $dataDir/DFree.out -time -node 3 4 -dof 1 2 3 disp;					# displacements of free nodes
recorder Node -file $dataDir/DBase.out -time -node 1 2 -dof 1 2 3 disp;					# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -node 1 2 -dof 1 2 3 reaction;					# support reaction
recorder Drift -file $dataDir/Drift.out -time -iNode 1 2 -jNode 3 4 -dof 1  -perpDirn 2 ;				# lateral drift
recorder Element -file $dataDir/FCol.out -time -ele 1 2 globalForce;						# element forces -- column
recorder Element -file $dataDir/FBeam.out -time -ele 3 globalForce;						# element forces -- beam
recorder Element -file $dataDir/ForceColSec1.out -time -ele 1 2 section 1 force;				# Column section forces, axial and moment, node i
recorder Element -file $dataDir/DefoColSec1.out -time -ele 1 2 section 1 deformation;				# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceColSec$numIntgrPts.out -time -ele 1 2 section $numIntgrPts force;		# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoColSec$numIntgrPts.out -time -ele 1 2 section $numIntgrPts deformation;		# section deformations, axial and curvature, node j
recorder Element -file $dataDir/ForceBeamSec1.out -time -ele 3 section 1 force;				# Beam section forces, axial and moment, node i
recorder Element -file $dataDir/DefoBeamSec1.out -time -ele 3 section 1 deformation;				# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceBeamSec$numIntgrPts.out -time -ele 3 section $numIntgrPts force;		# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoBeamSec$numIntgrPts.out -time -ele 3 section $numIntgrPts deformation;	# section deformations, axial and curvature, node j

# define GRAVITY -------------------------------------------------------------
set WzBeam [expr $Weight/$LBeam];
pattern Plain 1 Linear {
   eleLoad -ele 3 -type -beamUniform -$WzBeam ; # distributed superstructure-weight on beam
}
# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;			# convergence tolerance for test
constraints Plain;     			# how it handles boundary conditions
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
