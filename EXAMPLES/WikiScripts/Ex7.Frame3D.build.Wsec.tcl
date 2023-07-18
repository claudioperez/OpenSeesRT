# --------------------------------------------------------------------------------------------------
# Example 7. 3D Steel W-section Frame
#			Silvia Mazzoni & Frank McKenna, 2006
# nonlinearBeamColumn element, inelastic fiber section
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Output;			# set up name of data directory -- remove
file mkdir $dataDir; 			# create data directory
set GMdir "../GMfiles";		# ground-motion file directory
source LibUnits.tcl;			# define units
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel3D.tcl;		# procedure for displaying 3D perspectives of model
source WSection.tcl;		# procedure to define fiber W section

# define GEOMETRY -------------------------------------------------------------
# define structure-geometry paramters
set LCol [expr 12*$ft];		# column height (parallel to Y axis)
set LBeam [expr 20*$ft];		# beam length (parallel to X axis)
set LGird [expr 20*$ft];		# girder length (parallel to Z axis)

# ------ frame configuration
set NStory 3;			# number of stories above ground level
set NBay 1;			# number of bays in X direction
set NBayZ 1;			# number of bays in Z direction
puts "Number of Stories in Y: $NStory; Number of bays in X: $NBay; Number of bays in Z: $NBayZ"

# define NODAL COORDINATES
# calculate locations of beam/column intersections:
set X1 0.;
set X2 [expr $X1 + $LBeam];
set Y1 0.;
set Y2 [expr $Y1 + $LCol];
set Y3 [expr $Y2 + $LCol];
set Y4 [expr $Y3 + $LCol];
set Z1 0.0;
set Z2 [expr $Z1 + $LGird];

node 111 $X1 $Y1 $Z1;	# frame 1
node 112 $X2 $Y1 $Z1;
node 121 $X1 $Y2 $Z1;
node 122 $X2 $Y2 $Z1;
node 131 $X1 $Y3 $Z1;
node 132 $X2 $Y3 $Z1;
node 141 $X1 $Y4 $Z1;
node 142 $X2 $Y4 $Z1;
node 211 $X1 $Y1 $Z2;	# frame 2
node 212 $X2 $Y1 $Z2;
node 221 $X1 $Y2 $Z2;
node 222 $X2 $Y2 $Z2;
node 231 $X1 $Y3 $Z2;
node 232 $X2 $Y3 $Z2;
node 241 $X1 $Y4 $Z2;
node 242 $X2 $Y4 $Z2;


# define Rigid Floor Diaphragm
set RigidDiaphragm ON ;		# options: ON, OFF. specify this before the analysis parameters are set the constraints are handled differently.
set Xa [expr ($X2+$X1)/2];		# mid-span coordinate for rigid diaphragm
set Za [expr ($Z2+$Z1)/2];
# rigid-diaphragm nodes in center of each diaphram
node 1121 $Xa $Y2 $Za;		# master nodes for rigid diaphragm -- story 2, bay 1, frame 1-2
node 1131 $Xa $Y3 $Za;		# master nodes for rigid diaphragm -- story 3, bay 1, frame 1-2
node 1141 $Xa $Y4 $Za;		# master nodes for rigid diaphragm -- story 4, bay 1, frame 1-2
# Constraints for rigid diaphragm master nodes
fix 1121 0  1  0  1  0  1
fix 1131 0  1  0  1  0  1
fix 1141 0  1  0  1  0  1
# ------------------------define Rigid Diaphram, dof 2 is normal to floor
set perpDirn 2;
rigidDiaphragm $perpDirn 1121 121 122 221 222;	# level 2
rigidDiaphragm $perpDirn 1131 131 132 231 232;	# level 3 
rigidDiaphragm $perpDirn 1141 141 142 241 242;	# level 4 

# determine support nodes where ground motions are input, for multiple-support excitation
set iSupportNode "111 112 211 212"

# BOUNDARY CONDITIONS
fixY 0.0  1 1 1 0 1 0;		# pin all Y=0.0 nodes

# calculated MODEL PARAMETERS, particular to this model
# Set up parameters that are particular to the model for displacement control
set IDctrlNode 141;		# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement control
set LBuilding [expr $Y4];		# total building height

# Define  SECTIONS -------------------------------------------------------------
set SectionType FiberSection;		# options: Elastic FiberSection

# define section tags:
set ColSecTag 1
set BeamSecTag 2
set GirdSecTag 3
set ColSecTagFiber 4
set BeamSecTagFiber 5
set GirdSecTagFiber 6
set SecTagTorsion 70

if {$SectionType == "Elastic"} {
	# material properties:
	set Es [expr 29000*$ksi];		# Steel Young's Modulus
	set nu 0.3;			# Poisson's ratio
	set Gs [expr $Es/2./[expr 1+$nu]];  	# Torsional stiffness Modulus
	set J $Ubig;			# set large torsional stiffness
	# column sections: W27x114
	set AgCol [expr 33.5*pow($in,2)];		# cross-sectional area
	set IzCol [expr 4090.*pow($in,4)];		# moment of Inertia
	set IyCol [expr 159.*pow($in,4)];		# moment of Inertia
	# beam sections: W24x94
	set AgBeam [expr 27.7*pow($in,2)];		# cross-sectional area
	set IzBeam [expr 2700.*pow($in,4)];		# moment of Inertia
	set IyBeam [expr 109.*pow($in,4)];		# moment of Inertia
	# girder sections: W24x94
	set AgGird [expr 27.7*pow($in,2)];		# cross-sectional area
	set IzGird [expr 2700.*pow($in,4)];		# moment of Inertia
	set IyGird [expr 109.*pow($in,4)];		# moment of Inertia
	
	section Elastic $ColSecTag $Es $AgCol $IzCol $IyCol $Gs $J
	section Elastic $BeamSecTag $Es $AgBeam $IzBeam $IyBeam $Gs $J
	section Elastic $GirdSecTag $Es $AgGird $IzGird $IyGird $Gs $J

	set matIDhard 1;		# material numbers for recorder (this stressstrain recorder will be blank, as this is an elastic section)

} elseif {$SectionType == "FiberSection"} {
	# define MATERIAL properties 
	set Fy [expr 60.0*$ksi]
	set Es [expr 29000*$ksi];		# Steel Young's Modulus
	set nu 0.3;
	set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
	set Hiso 0
	set Hkin 1000
	set matIDhard 1
	uniaxialMaterial Hardening  $matIDhard $Es $Fy   $Hiso  $Hkin

	# ELEMENT properties
	# Structural-Steel W-section properties
	# column sections: W27x114
	set d [expr 27.29*$in];	# depth
	set bf [expr 10.07*$in];	# flange width
	set tf [expr 0.93*$in];	# flange thickness
	set tw [expr 0.57*$in];	# web thickness
	set nfdw 16;		# number of fibers along dw
	set nftw 2;		# number of fibers along tw
	set nfbf 16;		# number of fibers along bf
	set nftf 4;			# number of fibers along tf
	WSection  $ColSecTagFiber $matIDhard $d $bf $tf $tw $nfdw $nftw $nfbf $nftf
	# beam sections: W24x94
	set d [expr 24.31*$in];	# depth
	set bf [expr 9.065*$in];	# flange width
	set tf [expr 0.875*$in];	# flange thickness
	set tw [expr 0.515*$in];	# web thickness
	set nfdw 16;		# number of fibers along dw
	set nftw 2;		# number of fibers along tw
	set nfbf 16;		# number of fibers along bf
	set nftf 4;			# number of fibers along tf
	WSection  $BeamSecTagFiber $matIDhard $d $bf $tf $tw $nfdw $nftw $nfbf $nftf
	# girder sections: W24x94
	set d [expr 24.31*$in];	# depth
	set bf [expr 9.065*$in];	# flange width
	set tf [expr 0.875*$in];	# flange thickness
	set tw [expr 0.515*$in];	# web thickness
	set nfdw 16;		# number of fibers along dw
	set nftw 2;		# number of fibers along tw
	set nfbf 16;		# number of fibers along bf
	set nftf 4;			# number of fibers along tf
	WSection  $GirdSecTagFiber $matIDhard $d $bf $tf $tw $nfdw $nftw $nfbf $nftf
	
	# assign torsional Stiffness for 3D Model
	uniaxialMaterial Elastic $SecTagTorsion $Ubig
	section Aggregator $ColSecTag $SecTagTorsion T -section $ColSecTagFiber
	section Aggregator $BeamSecTag $SecTagTorsion T -section $BeamSecTagFiber
	section Aggregator $GirdSecTag $SecTagTorsion T -section $GirdSecTagFiber
} else {
	puts "No section has been defined"
	return -1
}

set QdlCol [expr 114*$lbf/$ft]; 		# W-section weight per length
set QBeam [expr 94*$lbf/$ft];		# W-section weight per length
set QGird [expr 94*$lbf/$ft];		# W-section weight per length

# define ELEMENTS -------------------------------------------------------
# set up geometric transformations of element
#   separate columns and beams, in case of P-Delta analysis for columns
#   in 3D model, assign vector vecxz
set IDColTransf 1; # all columns
set IDBeamTransf 2; # all beams
set IDGirdTransf 3; # all girders
set ColTransfType Linear ;				# options, Linear PDelta Corotational 
geomTransf $ColTransfType $IDColTransf  0 0 1 ; 	# only columns can have PDelta effects (gravity effects)
geomTransf Linear $IDBeamTransf 0 0 1
geomTransf Linear $IDGirdTransf 1 0 0

# Define Beam-Column Elements
set np 5;	# number of Gauss integration points for nonlinear curvature distribution

# Frame 1
# columns
element nonlinearBeamColumn 1111 111 121 $np $ColSecTag $IDColTransf;		# level 1-2
element nonlinearBeamColumn 1112 112 122 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 1121 121 131 $np $ColSecTag $IDColTransf;		# level 2-3
element nonlinearBeamColumn 1122 122 132 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 1131 131 141 $np $ColSecTag $IDColTransf;		# level 3-4
element nonlinearBeamColumn 1132 132 142 $np $ColSecTag $IDColTransf
# beams
element nonlinearBeamColumn 1221 121 122 $np $BeamSecTag $IDBeamTransf;	# level 2
element nonlinearBeamColumn 1231 131 132 $np $BeamSecTag $IDBeamTransf;	# level 3
element nonlinearBeamColumn 1241 141 142 $np $BeamSecTag $IDBeamTransf;	# level 4

# Frame 2
# columns
element nonlinearBeamColumn 2111 211 221 $np $ColSecTag $IDColTransf;		# level 1-2
element nonlinearBeamColumn 2112 212 222 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 2121 221 231 $np $ColSecTag $IDColTransf;		# level 2-3
element nonlinearBeamColumn 2122 222 232 $np $ColSecTag $IDColTransf
element nonlinearBeamColumn 2131 231 241 $np $ColSecTag $IDColTransf;		# level 3-4
element nonlinearBeamColumn 2132 232 242 $np $ColSecTag $IDColTransf
# beams
element nonlinearBeamColumn 2221 221 222 $np $BeamSecTag $IDBeamTransf;	# level 2
element nonlinearBeamColumn 2231 231 232 $np $BeamSecTag $IDBeamTransf;	# level 3
element nonlinearBeamColumn 2241 241 242 $np $BeamSecTag $IDBeamTransf;	# level 4

# girders connecting frames
# Frame 1-2
element nonlinearBeamColumn  1321 121 221 $np $GirdSecTag $IDGirdTransf;	# level 2
element nonlinearBeamColumn  1322 122 222 $np $GirdSecTag $IDGirdTransf;
element nonlinearBeamColumn  1331 131 231 $np $GirdSecTag $IDGirdTransf;	# level 3
element nonlinearBeamColumn  1332 132 232 $np $GirdSecTag $IDGirdTransf;
element nonlinearBeamColumn  1341 141 241 $np $GirdSecTag $IDGirdTransf;	# level 4
element nonlinearBeamColumn  1342 142 242 $np $GirdSecTag $IDGirdTransf;


# --------------------------------------------------------------------------------------------------------------------------------
# Define GRAVITY LOADS, weight and masses
# calculate dead load of frame, assume this to be an internal frame (do LL in a similar manner)
# calculate distributed weight along the beam length
set GammaConcrete [expr 150*$pcf];   			# Reinforced-Concrete weight density (weight per volume) 
set Tslab [expr 6*$in];			# 6-inch slab
set Lslab [expr $LGird/2]; 			# slab extends a distance of $LGird/2 in/out of plane
set DLfactor 1.0;				# scale dead load up a little
set Qslab [expr $GammaConcrete*$Tslab*$Lslab*$DLfactor]; 
set QdlBeam [expr $Qslab + $QBeam]; 	# dead load distributed along beam (one-way slab)
set QdlGird $QGird; 			# dead load distributed along girder
set WeightCol [expr $QdlCol*$LCol];  		# total Column weight
set WeightBeam [expr $QdlBeam*$LBeam]; 	# total Beam weight
set WeightGird [expr $QdlGird*$LGird]; 	# total Beam weight

# assign masses to the nodes that the columns are connected to 
# each connection takes the mass of 1/2 of each element framing into it (mass=weight/$g)
set Mmid  [expr ($WeightCol/2 + $WeightCol/2 +$WeightBeam/2+$WeightGird/2)/$g];
set Mtop  [expr ($WeightCol/2 + $WeightBeam/2+$WeightGird/2)/$g];

# frame 1
mass 121 $Mmid 0 $Mmid 0. 0. 0.;		# level 2
mass 122 $Mmid 0 $Mmid 0. 0. 0.;
mass 131 $Mmid 0 $Mmid 0. 0. 0.;		# level 3
mass 132 $Mmid 0 $Mmid 0. 0. 0.;
mass 141 $Mtop 0 $Mtop 0. 0. 0.;		# level 4
mass 142 $Mtop 0 $Mtop 0. 0. 0.;

# frame 2
mass 221 $Mmid 0 $Mmid 0. 0. 0.;		# level 2
mass 222 $Mmid 0 $Mmid 0. 0. 0.;
mass 231 $Mmid 0 $Mmid 0. 0. 0.;		# level 3
mass 232 $Mmid 0 $Mmid 0. 0. 0.;
mass 241 $Mtop 0 $Mtop 0. 0. 0.;		# level 4
mass 242 $Mtop 0 $Mtop 0. 0. 0.;

set FloorWeight2 [expr 4*$WeightCol + 2*$WeightGird + 2*$WeightBeam]
set FloorWeight3 [expr 4*$WeightCol + 2*$WeightGird + 2*$WeightBeam]
set FloorWeight4 [expr 2*$WeightCol + 2*$WeightGird + 2*$WeightBeam]
set WeightTotal [expr $FloorWeight2+$FloorWeight3+$FloorWeight4];			# total building weight
set MassTotal [expr $WeightTotal/$g];							# total building mass

# --------------------------------------------------------------------------------------------------------------------------------
# LATERAL-LOAD distribution for static pushover analysis
# calculate distribution of lateral load based on mass/weight distributions along building height
# Fj = WjHj/sum(WiHi)  * Weight   at each floor j
set sumWiHi [expr $FloorWeight2*$Y2+$FloorWeight3*$Y3+$FloorWeight4*$Y4]; 		# sum of storey weight times height, for lateral-load distribution
set WiHi2 [expr $FloorWeight2*$Y2]; 		# storey weight times height, for lateral-load distribution
set WiHi3 [expr $FloorWeight3*$Y3]; 		# storey weight times height, for lateral-load distribution
set WiHi4 [expr $FloorWeight4*$Y4]; 		# storey weight times height, for lateral-load distribution
set F2 [expr $WiHi2/$sumWiHi*$WeightTotal];	# lateral load at level
set F3 [expr $WiHi3/$sumWiHi*$WeightTotal];	# lateral load at level
set F4 [expr $WiHi4/$sumWiHi*$WeightTotal];	# lateral load at level


# Define RECORDERS -------------------------------------------------------------
recorder Node -file $dataDir/DFree.out -time -node 141 -dof 1 2 3 disp;			# displacements of free node
recorder Node -file $dataDir/DBase.out -time -node 111 112 211 212   -dof 1 2 3 disp;		# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -node 111 112 211 212   -dof 1 2 3 reaction;		# support reaction
recorder Drift -file $dataDir/DrNode.out -time -iNode 111 -jNode 141 -dof 1 -perpDirn 2;		# lateral drift
recorder Element -file $dataDir/Fel1.out -time -ele 1111 localForce;				# element forces in local coordinates
recorder Element -xml $dataDir/PlasticRotation1.out -time -ele 1111 plasticRotation;				# element forces in local coordinates
recorder Element -file $dataDir/ForceEle1sec1.out -time -ele 1111 section 1 force;			# section forces, axial and moment, node i
recorder Element -file $dataDir/DefoEle1sec1.out -time -ele 11111 section 1 deformation;		# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceEle1sec$np.out -time -ele 111 section $np force;			# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoEle1sec$np.out -time -ele 1111 section $np deformation;		# section deformations, axial and curvature, node j
set yFiber [expr 0.];								# fiber location for stress-strain recorder, local coords
set zFiber [expr 0.];								# fiber location for stress-strain recorder, local coords
recorder Element -file $dataDir/SSreinfEle1sec1.out -time -ele 1111 section $np fiber $yFiber $zFiber stressStrain;	# steel fiber stress-strain, node i

# Define DISPLAY -------------------------------------------------------------
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set dAmp 2;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel3D NodeNumbers $dAmp $xLoc1 $yLoc1  $xPixels $yPixels

# define GRAVITY -------------------------------------------------------------
# GRAVITY LOADS # define gravity load applied to beams and columns -- 	eleLoad applies loads in local coordinate axis
pattern Plain 101 Linear {
# Frame 1
# columns
	eleLoad -ele 1111 -type -beamUniform 0. 0. -$QdlCol;		# level 1-2
	eleLoad -ele 1112 -type -beamUniform 0. 0. -$QdlCol
	eleLoad -ele 1121 -type -beamUniform 0. 0. -$QdlCol;		# level 2-3
	eleLoad -ele 1122 -type -beamUniform 0. 0. -$QdlCol
	eleLoad -ele 1131 -type -beamUniform 0. 0. -$QdlCol;		# level 3-4
	eleLoad -ele 1132 -type -beamUniform 0. 0. -$QdlCol
# beams
	eleLoad -ele 1221 -type -beamUniform -$QdlBeam 0.;		# level 2
	eleLoad -ele 1231 -type -beamUniform -$QdlBeam 0.;		# level 3
	eleLoad -ele 1241 -type -beamUniform -$QdlBeam 0.;		# level 4

# Frame 2
# columns
	eleLoad -ele 2111 -type -beamUniform 0. 0. -$QdlCol;		# level 1-2
	eleLoad -ele 2112 -type -beamUniform 0. 0. -$QdlCol
	eleLoad -ele 2121 -type -beamUniform 0. 0. -$QdlCol;		# level 2-3
	eleLoad -ele 2122 -type -beamUniform 0. 0. -$QdlCol
	eleLoad -ele 2131 -type -beamUniform 0. 0. -$QdlCol;		# level 3-4
	eleLoad -ele 2132 -type -beamUniform 0. 0. -$QdlCol
# beams
	eleLoad -ele 2221 -type -beamUniform -$QdlBeam 0.;		# level 2
	eleLoad -ele 2231 -type -beamUniform -$QdlBeam 0.;		# level 3
	eleLoad -ele 2241 -type -beamUniform -$QdlBeam 0.;		# level 4

# girders connecting frames
# Frame 1-2
	eleLoad -ele 1321 -type -beamUniform -$QdlGird 0.;		# level 2
	eleLoad -ele 1322 -type -beamUniform -$QdlGird 0.;
	eleLoad -ele 1331 -type -beamUniform -$QdlGird 0.;		# level 3
	eleLoad -ele 1332 -type -beamUniform -$QdlGird 0.;
	eleLoad -ele 1341 -type -beamUniform -$QdlGird 0.;		# level 4
	eleLoad -ele 1342 -type -beamUniform -$QdlGird 0.;
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
test EnergyIncr $Tol 6 ; 		# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;			# use Newton's solution algorithm: updates tangent stiffness at every iteration
set NstepGravity 10;  		# apply gravity in 10 steps
set DGravity [expr 1./$NstepGravity]; 	# first load increment;
integrator LoadControl $DGravity;	# determine the next time step for an analysis
analysis Static;			# define type of analysis static or transient
analyze $NstepGravity;		# apply gravity

# ------------------------------------------------- maintain constant gravity loads and reset time to zero
loadConst -time 0.0
set Tol 1.0e-6;			# reduce tolerance after gravity loads
puts "Model Built"


