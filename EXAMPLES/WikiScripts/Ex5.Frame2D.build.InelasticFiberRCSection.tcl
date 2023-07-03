# --------------------------------------------------------------------------------------------------
# Example 5. 2D Frame --  Build Model
# nonlinearBeamColumn element, inelastic fiber section -- Reinforced Concrete Section
#			Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory (you can remove this)
file mkdir $dataDir; 			# create data directory
set GMdir "../GMfiles/";		# ground-motion file directory
source LibUnits.tcl;			# define units
source BuildRCrectSection.tcl;		# procedure for definining RC fiber section

# define GEOMETRY -------------------------------------------------------------
# define structure-geometry paramters
set LCol [expr 14*$ft];	# column height
set LBeam [expr 24*$ft];	# beam length

# calculate locations of beam/column intersections:
set X1 0.;
set X2 [expr $X1 + $LBeam];
set X3 [expr $X2 + $LBeam];
set X4 [expr $X3 + $LBeam];
set Y1 0.;
set Y2 [expr $Y1 + $LCol];
set Y3 [expr $Y2 + $LCol];
set Y4 [expr $Y3 + $LCol];

# define nodal coordinates
node 11 $X1 $Y1
node 12 $X2 $Y1
node 13 $X3 $Y1
node 14 $X4 $Y1
node 21 $X1 $Y2
node 22 $X2 $Y2
node 23 $X3 $Y2
node 24 $X4 $Y2
node 31 $X1 $Y3
node 32 $X2 $Y3
node 33 $X3 $Y3
node 34 $X4 $Y3
node 41 $X1 $Y4
node 42 $X2 $Y4
node 43 $X3 $Y4
node 44 $X4 $Y4

# Set up parameters that are particular to the model for displacement control
set IDctrlNode 41;		# node where displacement is read for displacement control
set IDctrlDOF 1;		# degree of freedom of displacement read for displacement control
set NStory 3;		# number of stories above ground level
set NBay 3;		# number of bays
set LBuilding $Y4;		# total building height

# BOUNDARY CONDITIONS
fix 11 1 1 0
fix 12 1 1 0
fix 13 1 1 0
fix 14 1 1 0

# Define SECTIONS -------------------------------------------------------------
set SectionType FiberSection;		# options: Elastic FiberSection

# define section tags:
set ColSecTag 1
set BeamSecTag 2

# Section Properties:
set HCol [expr 24*$in];		# square-Column width
set BCol $HCol
set HBeam [expr 42*$in];		# Beam depth -- perpendicular to bending axis
set BBeam [expr 24*$in];		# Beam width -- parallel to bending axis

if {$SectionType == "Elastic"} {
	# material properties:
	set fc 4000*$psi;			# concrete nominal compressive strength
	set Ec [expr 57*$ksi*pow($fc/$psi,0.5)];	# concrete Young's Modulus
	# column section properties:
	set AgCol [expr $HCol*$BCol];		# rectuangular-Column cross-sectional area
	set IzCol [expr 0.5*1./12*$BCol*pow($HCol,3)];	# about-local-z Rect-Column gross moment of inertial
	# beam sections:
	set AgBeam [expr $HBeam*$BBeam];		# rectuangular-Beam cross-sectional area
	set IzBeam [expr 0.5*1./12*$BBeam*pow($HBeam,3)];	# about-local-z Rect-Beam cracked moment of inertial
		
	section Elastic $ColSecTag $Ec $AgCol $IzCol 
	section Elastic $BeamSecTag $Ec $AgBeam $IzBeam 

} elseif {$SectionType == "FiberSection"} {
	# MATERIAL parameters 
	source LibMaterialsRC.tcl;	# define library of Reinforced-concrete Materials

	# FIBER SECTION properties 
	# Column section geometry:
	set cover [expr 2.5*$in];	# rectangular-RC-Column cover
	set numBarsTopCol 8;		# number of longitudinal-reinforcement bars on top layer
	set numBarsBotCol 8;		# number of longitudinal-reinforcement bars on bottom layer
	set numBarsIntCol 6;		# TOTAL number of reinforcing bars on the intermediate layers
	set barAreaTopCol [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaBotCol [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaIntCol [expr 1.*$in*$in];	# longitudinal-reinforcement bar area

	set numBarsTopBeam 6;		# number of longitudinal-reinforcement bars on top layer
	set numBarsBotBeam 6;		# number of longitudinal-reinforcement bars on bottom layer
	set numBarsIntBeam 2;		# TOTAL number of reinforcing bars on the intermediate layers
	set barAreaTopBeam [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaBotBeam [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaIntBeam [expr 1.*$in*$in];	# longitudinal-reinforcement bar area

	set nfCoreY 20;		# number of fibers in the core patch in the y direction
	set nfCoreZ 20;		# number of fibers in the core patch in the z direction
	set nfCoverY 20;		# number of fibers in the cover patches with long sides in the y direction
	set nfCoverZ 20;		# number of fibers in the cover patches with long sides in the z direction
	# rectangular section with one layer of steel evenly distributed around the perimeter and a confined core.
	BuildRCrectSection $ColSecTag $HCol $BCol $cover $cover $IDconcCore  $IDconcCover $IDSteel $numBarsTopCol $barAreaTopCol $numBarsBotCol $barAreaBotCol $numBarsIntCol $barAreaIntCol  $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ
	BuildRCrectSection $BeamSecTag $HBeam $BBeam $cover $cover $IDconcCore  $IDconcCover $IDSteel $numBarsTopBeam $barAreaTopBeam $numBarsBotBeam $barAreaBotBeam $numBarsIntBeam $barAreaIntBeam  $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ
} else {
	puts "No section has been defined"
	return -1
}


# define ELEMENTS
# set up geometric transformations of element
#   separate columns and beams, in case of P-Delta analysis for columns
set IDColTransf 1; # all columns
set IDBeamTransf 2; # all beams
set ColTransfType Linear ;			# options, Linear PDelta Corotational 
geomTransf $ColTransfType $IDColTransf  ; 	# only columns can have PDelta effects (gravity effects)
geomTransf Linear $IDBeamTransf

# Define Beam-Column Elements
set np 5;	# number of Gauss integration points for nonlinear curvature distribution-- np=2 for linear distribution ok
# columns
element nonlinearBeamColumn 111 11 21 $np $ColSecTag $IDColTransf;		# level 1-2
element nonlinearBeamColumn 112 12 22 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 113 13 23 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 114 14 24 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 121 21 31 $np $ColSecTag $IDColTransf;		# level 2-3
element nonlinearBeamColumn 122 22 32 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 123 23 33 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 124 24 34 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 131 31 41 $np $ColSecTag $IDColTransf;		# level 3-4
element nonlinearBeamColumn 132 32 42 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 133 33 43 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 134 34 44 $np $ColSecTag $IDColTransf
# beams
element nonlinearBeamColumn 221 21 22 $np $BeamSecTag $IDBeamTransf;		# level 2
element nonlinearBeamColumn 222 22 23 $np $BeamSecTag $IDBeamTransf;
element nonlinearBeamColumn 223 23 24 $np $BeamSecTag $IDBeamTransf;
element nonlinearBeamColumn 231 31 32 $np $BeamSecTag $IDBeamTransf;		# level 3
element nonlinearBeamColumn 232 32 33 $np $BeamSecTag $IDBeamTransf;
element nonlinearBeamColumn 233 33 34 $np $BeamSecTag $IDBeamTransf;
element nonlinearBeamColumn 241 41 42 $np $BeamSecTag $IDBeamTransf;		# level 4
element nonlinearBeamColumn 242 42 43 $np $BeamSecTag $IDBeamTransf;
element nonlinearBeamColumn 243 43 44 $np $BeamSecTag $IDBeamTransf;
	
# Define GRAVITY LOADS, weight and masses
# calculate dead load of frame, assume this to be an internal frame (do LL in a similar manner)
# calculate distributed weight along the beam length
set GammaConcrete [expr 150*$pcf];   		# Reinforced-Concrete floor slabs
set Tslab [expr 6*$in];			# 6-inch slab
set Lslab [expr 2*$LBeam/2]; 			# assume slab extends a distance of $LBeam1/2 in/out of plane
set Qslab [expr $GammaConcrete*$Tslab*$Lslab]; 
set QBeam [expr 94*$lbf/$ft];		# W-section weight per length
set QdlBeam [expr $Qslab + $QBeam]; 	# dead load distributed along beam.
set QdlCol [expr 114*$lbf/$ft]; 	# W-section weight per length
set WeightCol [expr $QdlCol*$LCol];  		# total Column weight
set WeightBeam [expr $QdlBeam*$LBeam]; 	# total Beam weight

# assign masses to the nodes that the columns are connected to 
# each connection takes the mass of 1/2 of each element framing into it (mass=weight/$g)
mass 21   [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;			# level 2
mass 22  [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 23 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 24 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;
mass 31 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;			# level 3
mass 32 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 33 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 34 [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;
mass 41 [expr ($WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;					# level 4
mass 42 [expr ($WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 43 [expr ($WeightCol/2 +$WeightBeam/2 +$WeightBeam/2)/$g] 0. 0.;
mass 44 [expr ($WeightCol/2 +$WeightBeam/2)/$g] 0. 0.;
# calculate total Floor Mass
set WeightFloor2 [expr $WeightCol*4/2+$WeightCol*4/2+3*$WeightBeam];			# level 2 weight
set WeightFloor3 [expr $WeightCol*4/2+$WeightCol*4/2+3*$WeightBeam];
set WeightFloor4 [expr $WeightCol*4/2+3*$WeightBeam];
set WeightTotal [expr $WeightFloor2 + $WeightFloor3 + $WeightFloor4];			# total frame weight
set MassFloor2 [expr $WeightFloor2/$g];
set MassFloor3 [expr $WeightFloor3/$g];
set MassFloor4 [expr $WeightFloor4/$g];
set MassTotal [expr $MassFloor2+$MassFloor3+$MassFloor4];				# total frame mass

# LATERAL-LOAD distribution for static pushover analysis
# calculate distribution of lateral load based on mass/weight distributions along building height
# Fj = WjHj/sum(WiHi)  * Weight   at each floor j
set sumWiHi [expr $WeightFloor2*$Y2 + $WeightFloor3*$Y3 + $WeightFloor4*$Y4];	# denominator
set Fj2 [expr $WeightFloor2*$Y2/$sumWiHi*$WeightTotal];		# total for floor 2
set Fj3 [expr $WeightFloor3*$Y3/$sumWiHi*$WeightTotal];		# total for floor 3
set Fj4 [expr $WeightFloor4*$Y4/$sumWiHi*$WeightTotal];		# total for floor 4
set Fi2 [expr $Fj2/4];			# per node on floor 2
set Fi3 [expr $Fj3/4];			# per node on floor 3
set Fi4 [expr $Fj4/4];			# per node on floor 4
set iFi "$Fi2 $Fi3 $Fi4";			# vectorize

# Define RECORDERS -------------------------------------------------------------
recorder Node -file $dataDir/DFree.out -time -node 41 -dof 1 2 3 disp;			# displacements of free node
recorder Node -file $dataDir/DBase.out -time -node 11 12 13 14 -dof 1 2 3 disp;		# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -node 11 12 13 14 -dof 1 2 3 reaction;		# support reaction
recorder Drift -file $dataDir/DrNode.out -time -iNode 11 -jNode 41 -dof 1 -perpDirn 2;		# lateral drift
recorder Element -file $dataDir/Fel1.out -time -ele 111 localForce;				# element forces in local coordinates
recorder Element -file $dataDir/ForceEle1sec1.out -time -ele 111 section 1 force;			# section forces, axial and moment, node i
recorder Element -file $dataDir/DefoEle1sec1.out -time -ele 111 section 1 deformation;			# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceEle1sec$np.out -time -ele 111 section $np force;			# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoEle1sec$np.out -time -ele 111 section $np deformation;		# section deformations, axial and curvature, node j
recorder Element -file $dataDir/SSEle1sec1.out -time -ele 111 section $np fiber 0 0 $IDSteel  stressStrain;	# steel fiber stress-strain, node i

# define GRAVITY -------------------------------------------------------------
# GRAVITY LOADS # define gravity load applied to beams and columns -- eleLoad applies loads in local coordinate axis
pattern Plain 101 Linear {
   eleLoad -ele 221 222 223  -type -beamUniform -$QdlBeam; ; # beams level 2 (in -ydirection)
   eleLoad -ele 231 232 233 -type -beamUniform -$QdlBeam; 
   eleLoad -ele 241 242 243  -type -beamUniform -$QdlBeam
   eleLoad -ele 111 112 113 114 -type -beamUniform 0 -$QdlCol; # columns level 1-2 (in -xdirection)
   eleLoad -ele 121 122 123 124 -type -beamUniform 0 -$QdlCol; 
   eleLoad -ele 131 132 133 134 -type -beamUniform 0 -$QdlCol; 
}
# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;			# convergence tolerance for test
variable constraintsTypeGravity Plain;		# default;
if {  [info exists RigidDiaphragm] == 1} {
	if {$RigidDiaphragm=="ON"} {
		variable constraintsTypeGravity Lagrange;	#  large model: try Transformation
	};	# if rigid diaphragm is on
};	# if rigid diaphragm exists
constraints $constraintsTypeGravity ;     		# how it handles boundary conditions
numberer RCM;			# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral ;		# how to store and solve the system of equations in the analysis (large model: try UmfPack)
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


