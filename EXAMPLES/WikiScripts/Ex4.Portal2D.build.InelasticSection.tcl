# --------------------------------------------------------------------------------------------------
# Example4. 2D Portal Frame--  Build Model
# nonlinearBeamColumn element, uniaxial inelastic section
#			Silvia Mazzoni & Frank McKenna, 2006
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
set GMdir "GMfiles";			# ground-motion file directory
source LibUnits.tcl;			# define basic and system units

# define GEOMETRY -------------------------------------------------------------
set LCol [expr 36*$ft]; 		# column length
set LBeam [expr 42*$ft];		# beam length
set Weight [expr 4000.*$kip]; 		# superstructure weight
# define section geometry
set HCol [expr 5.*$ft]; 		# Column Depth
set BCol [expr 4.*$ft];		# Column Width
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

# Define ELEMENTS & SECTIONS  -------------------------------------------------------------
set ColMatTagFlex 2;			# assign a tag number to the column flexural behavior
set ColMatTagAxial 3;			# assign a tag number to the column axial behavior	
set ColSecTag 1;				# assign a tag number to the column section tag
set BeamSecTag 2;				# assign a tag number to the beam section tag

# MATERIAL parameters
set fc [expr -4*$ksi]; 		# CONCRETE Compressive Strength (+Tension, -Compression)
set Ec [expr 57*$ksi*sqrt(-$fc/$psi)]; 	# Concrete Elastic Modulus

# COLUMN section
# calculated stiffness parameters
set EICol [expr $Ec*$IzCol];			# EI, for moment-curvature relationship
set EACol [expr $Ec*$ACol];			# EA, for axial-force-strain relationship
set MyCol [expr 130000*$kip*$in];		# yield moment
set PhiYCol [expr 0.65e-4/$in];			# yield curvature
set EIColCrack [expr $MyCol/$PhiYCol];	# cracked section inertia
set b 0.01 ;				# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $ColMatTagFlex $MyCol $EIColCrack $b; 		# bilinear behavior for flexure
uniaxialMaterial Elastic $ColMatTagAxial $EACol;				# this is not used as a material, this is an axial-force-strain response
section Aggregator $ColSecTag $ColMatTagAxial P $ColMatTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)

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
