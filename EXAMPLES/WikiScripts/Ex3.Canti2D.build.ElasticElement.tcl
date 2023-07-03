# --------------------------------------------------------------------------------------------------
# Example 3. 2D Cantilever -- Build Model
# elasticBeamColumn element
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
wipe;					# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;		# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;				# set up name for data directory
file mkdir $dataDir/; 			# create data directory
set GMdir "../GMfiles";			# ground-motion file directory

# define UNITS ----------------------------------------------------------------------------
set in 1.; 				# define basic units -- output units
set kip 1.; 			# define basic units -- output units
set sec 1.; 			# define basic units -- output units
set LunitTXT "inch";			# define basic-unit text for output
set FunitTXT "kip";			# define basic-unit text for output
set TunitTXT "sec";			# define basic-unit text for output
set ft [expr 12.*$in]; 		# define engineering units
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set lbf [expr $psi*$in*$in];		# pounds force
set pcf [expr $lbf/pow($ft,3)];		# pounds per cubic foot
set in2 [expr $in*$in]; 		# inch^2
set in4 [expr $in*$in*$in*$in]; 		# inch^4
set cm [expr $in/2.54];		# centimeter, needed for displacement input in MultipleSupport excitation
set PI [expr 2*asin(1.0)]; 		# define constants
set g [expr 32.2*$ft/pow($sec,2)]; 	# gravitational acceleration
set Ubig 1.e10; 			# a really large number
set Usmall [expr 1/$Ubig]; 		# a really small number

# define GEOMETRY -------------------------------------------------------------
set LCol [expr 36*$ft]; 		# column length
set Weight [expr 2000.*$kip]; 		# superstructure weight
# define section geometry
set HCol [expr 5.*$ft]; 		# Column Depth
set BCol [expr 5.*$ft];		# Column Width

# calculated parameters
set PCol [expr $Weight]; 		# nodal dead-load weight per column
set Mass [expr $PCol/$g];		# nodal mass
# calculated geometry parameters
set ACol [expr $BCol*$HCol];					# cross-sectional area
set IzCol [expr 1./12.*$BCol*pow($HCol,3)]; 			# Column moment of inertia

# nodal coordinates:
node 1 0 0;			# node#, X, Y
node 2 0 $LCol 		

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			# node DX DY RZ

# we need to set up parameters that are particular to the model.
set IDctrlNode 2;			# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement control
set iSupportNode "1";		# define support node, if needed.

# nodal masses:
mass 2 $Mass  1e-9 0.;		# node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

# Define ELEMENTS -------------------------------------------------------------
# Material parameters
set fc [expr -4.*$ksi]; 		# CONCRETE Compressive Strength (+Tension, -Compression)
set Ec [expr 57*$ksi*sqrt(-$fc/$psi)]; 	# Concrete Elastic Modulus

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set ColTransfTag 1; 			# associate a tag to column transformation
set ColTransfType Linear ;			# options, Linear PDelta Corotational 
geomTransf $ColTransfType $ColTransfTag ; 	

# element connectivity:
element elasticBeamColumn 1 1 2 $ACol $Ec $IzCol $ColTransfTag;			# self-explanatory when using variables

# Define RECORDERS -------------------------------------------------------------
recorder Node -file $dataDir/DFree.out -time -node 2 -dof 1 2 3 disp;		# displacements of free nodes
recorder Node -file $dataDir/DBase.out -time -node 1 -dof 1 2 3 disp;		# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -node 1 -dof 1 2 3 reaction;		# support reaction
recorder Drift -file $dataDir/Drift.out -time -iNode 1 -jNode 2 -dof 1   -perpDirn 2 ;	# lateral drift
recorder Element -file $dataDir/FCol.out -time -ele 1 globalForce;			# element forces -- column
recorder Element -xml $dataDir/PlasticRotation.out -time -ele 1 plasticRotation;		# section deformations, axial and curvature, node j

# define GRAVITY -------------------------------------------------------------
pattern Plain 1 Linear {
   load 2 0 -$PCol 0
}

# ------------------------------------------------- apply gravity load
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
