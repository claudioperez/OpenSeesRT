# --------------------------------------------------------------------------------------------------
# Example 8. 3D RC frame
#		Silvia Mazzoni & Frank McKenna, 2006
# nonlinearBeamColumn element, inelastic fiber section
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Output;			# set up name of data directory
file mkdir $dataDir; 			# create data directory
set GMdir "GMfiles";		# ground-motion file directory
source LibUnits_IR.tcl;			# define units
source DisplayPlane_IR.tcl;		# procedure for displaying a plane in model
source DisplayModel3D_IR.tcl;		# procedure for displaying 3D perspectives of model
source Library/BuildRCrectSection_IR.tcl;		# procedure for definining RC fiber section

# ------ frame configuration
set NStory 2;			# number of stories above ground level
set NBay 2;			# number of bays in X direction
set NBayZ 3;			# number of bays in Z direction
puts "Number of Stories in Y: $NStory Number of bays in X: $NBay Number of bays in Z: $NBayZ"
set NFrame [expr $NBayZ + 1];	# actually deal with frames in Z direction, as this is an easy extension of the 2d model

# define GEOMETRY -------------------------------------------------------------
# define structure-geometry paramters
set LCol [expr 14*$ft];		# column height (parallel to Y axis)
set LBeam [expr 24*$ft];		# beam length (parallel to X axis)
set LGird [expr 24*$ft];		# girder length (parallel to Z axis)

# define NODAL COORDINATES
set Dlevel 10000;	# numbering increment for new-level nodes
set Dframe 100;	# numbering increment for new-frame nodes
for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
	set Z [expr ($frame-1)*$LGird];
	for {set level 1} {$level <=[expr $NStory+1]} {incr level 1} {
		set Y [expr ($level-1)*$LCol];
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
			set X [expr ($pier-1)*$LBeam];
			set nodeID [expr $level*$Dlevel+$frame*$Dframe+$pier]
			node $nodeID $X $Y $Z;		# actually define node
		}
	}
}

# rigid diaphragm nodes
set RigidDiaphragm ON ;		# options: ON, OFF. specify this before the analysis parameters are set the constraints are handled differently.
set Xa [expr ($NBay*$LBeam)/2];		# mid-span coordinate for rigid diaphragm
set Za [expr ($NFrame-1)*$LGird/2];
set iMasterNode ""
for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
	set Y [expr ($level-1)*$LCol];
	# rigid-diaphragm nodes in center of each diaphram
	set MasterNodeID [expr 9900+$level]
	node $MasterNodeID $Xa $Y $Za;		# master nodes for rigid diaphragm
	fix $MasterNodeID 0  1  0  1  0  1;		# constrain other dofs that don't belong to rigid diaphragm control
	lappend iMasterNode $MasterNodeID
	set perpDirn 2;				# perpendicular to plane of rigid diaphragm
	for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
			set nodeID [expr $level*$Dlevel+$frame*$Dframe+$pier]
			rigidDiaphragm $perpDirn $MasterNodeID $nodeID; 	# define Rigid Diaphram,
		}
	}
}

# determine support nodes where ground motions are input, for multiple-support excitation
set iSupportNode ""
for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
	set level 1
	for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
		set nodeID [expr $level*$Dlevel+$frame*$Dframe+$pier]
		lappend iSupportNode $nodeID
	}
}

# BOUNDARY CONDITIONS
fixY 0.0  1 1 1 0 1 0;		# pin all Y=0.0 nodes

# calculated MODEL PARAMETERS, particular to this model
# Set up parameters that are particular to the model for displacement control
set IDctrlNode [expr int(($NStory+1)*$Dlevel+(1*$Dframe)+1)];		# node where displacement is read for displacement control
set IDctrlDOF 1;					# degree of freedom of displacement read for displacement control
set LBuilding [expr $NStory*$LCol];			# total building height

# Define SECTIONS -------------------------------------------------------------
set SectionType FiberSection;		# options: Elastic FiberSection

# define section tags:
set ColSecTag 1
set BeamSecTag 2
set GirdSecTag 3
set ColSecTagFiber 4
set BeamSecTagFiber 5
set GirdSecTagFiber 6
set SecTagTorsion 70

# Section Properties:
set HCol [expr 18*$in];		# square-Column width
set BCol $HCol
set HBeam [expr 24*$in];		# Beam depth -- perpendicular to bending axis
set BBeam [expr 18*$in];		# Beam width -- parallel to bending axis
set HGird [expr 24*$in];		# Girder depth -- perpendicular to bending axis
set BGird [expr 18*$in];		# Girder width -- parallel to bending axis

if {$SectionType == "Elastic"} {
	# material properties:
	set fc 4000*$psi;			# concrete nominal compressive strength
	set Ec [expr 57*$ksi*pow($fc/$psi,0.5)];	# concrete Young's Modulus
	set nu 0.2;			# Poisson's ratio
	set Gc [expr $Ec/2./[expr 1+$nu]];  	# Torsional stiffness Modulus
	set J $Ubig;			# set large torsional stiffness
	# column section properties:
	set AgCol [expr $HCol*$BCol];		# rectuangular-Column cross-sectional area
	set IzCol [expr 0.5*1./12*$BCol*pow($HCol,3)];	# about-local-z Rect-Column gross moment of inertial
	set IyCol [expr 0.5*1./12*$HCol*pow($BCol,3)];	# about-local-z Rect-Column gross moment of inertial
	# beam sections:
	set AgBeam [expr $HBeam*$BBeam];		# rectuangular-Beam cross-sectional area
	set IzBeam [expr 0.5*1./12*$BBeam*pow($HBeam,3)];	# about-local-z Rect-Beam cracked moment of inertial
	set IyBeam [expr 0.5*1./12*$HBeam*pow($BBeam,3)];	# about-local-y Rect-Beam cracked moment of inertial
	# girder sections:
	set AgGird [expr $HGird*$BGird];		# rectuangular-Girder cross-sectional area
	set IzGird [expr 0.5*1./12*$BGird*pow($HGird,3)];	# about-local-z Rect-Girder cracked moment of inertial
	set IyGird [expr 0.5*1./12*$HGird*pow($BGird,3)];	# about-local-y Rect-Girder cracked moment of inertial
		
	section Elastic $ColSecTag $Ec $AgCol $IzCol $IyCol $Gc $J
	section Elastic $BeamSecTag $Ec $AgBeam $IzBeam $IyBeam $Gc $J
	section Elastic $GirdSecTag $Ec $AgGird $IzGird $IyGird $Gc $J

	set IDconcCore  1;		# material numbers for recorder (this stressstrain recorder will be blank, as this is an elastic section)
	set IDSteel  2;			# material numbers for recorder (this stressstrain recorder will be blank, as this is an elastic section)

} elseif {$SectionType == "FiberSection"} {
	# MATERIAL parameters 
	source LibMaterialsRC_IR.tcl;	# define library of Reinforced-concrete Materials
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

	set numBarsTopGird 6;		# number of longitudinal-reinforcement bars on top layer
	set numBarsBotGird 6;		# number of longitudinal-reinforcement bars on bottom layer
	set numBarsIntGird 2;		# TOTAL number of reinforcing bars on the intermediate layers
	set barAreaTopGird [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaBotGird [expr 1.*$in*$in];	# longitudinal-reinforcement bar area
	set barAreaIntGird [expr 1.*$in*$in];	# longitudinal-reinforcement bar area

	set nfCoreY 12;		# number of fibers in the core patch in the y direction
	set nfCoreZ 12;		# number of fibers in the core patch in the z direction
	set nfCoverY 8;		# number of fibers in the cover patches with long sides in the y direction
	set nfCoverZ 8;		# number of fibers in the cover patches with long sides in the z direction
	# rectangular section with one layer of steel evenly distributed around the perimeter and a confined core.
	BuildRCrectSection $ColSecTagFiber $HCol $BCol $cover $cover $IDconcCore  $IDconcCover $IDSteel $numBarsTopCol $barAreaTopCol $numBarsBotCol $barAreaBotCol $numBarsIntCol $barAreaIntCol  $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ
	BuildRCrectSection $BeamSecTagFiber $HBeam $BBeam $cover $cover $IDconcCore  $IDconcCover $IDSteel $numBarsTopBeam $barAreaTopBeam $numBarsBotBeam $barAreaBotBeam $numBarsIntBeam $barAreaIntBeam  $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ
	BuildRCrectSection $GirdSecTagFiber $HGird $BGird $cover $cover $IDconcCore  $IDconcCover $IDSteel $numBarsTopGird $barAreaTopGird $numBarsBotGird $barAreaBotGird $numBarsIntGird $barAreaIntGird  $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ

	# assign torsional Stiffness for 3D Model
	uniaxialMaterial Elastic $SecTagTorsion $Ubig
	section Aggregator $ColSecTag $SecTagTorsion T -section $ColSecTagFiber
	section Aggregator $BeamSecTag $SecTagTorsion T -section $BeamSecTagFiber
	section Aggregator $GirdSecTag $SecTagTorsion T -section $GirdSecTagFiber
} else {
	puts "No section has been defined"
	return -1
}
set GammaConcrete [expr 150*$pcf];
set QdlCol [expr $GammaConcrete*$HCol*$BCol];	# self weight of Column, weight per length
set QBeam [expr $GammaConcrete*$HBeam*$BBeam];	# self weight of Beam, weight per length
set QGird [expr $GammaConcrete*$HGird*$BGird];	# self weight of Gird, weight per length

# --------------------------------------------------------------------------------------------------------------------------------
# define ELEMENTS
# set up geometric transformations of element
#   separate columns and beams, in case of P-Delta analysis for columns
set IDColTransf 1; # all columns
set IDBeamTransf 2; # all beams
set IDGirdTransf 3; # all girds
set ColTransfType Linear ;		# options for columns: Linear PDelta  Corotational 
geomTransf $ColTransfType  $IDColTransf  0 0 1;			# orientation of column stiffness affects bidirectional response.
geomTransf Linear $IDBeamTransf 0 0 1
geomTransf Linear $IDGirdTransf 1 0 0

# Define Beam-Column Elements
set numIntgrPts 5;	# number of Gauss integration points for nonlinear curvature distribution
# columns
set N0col [expr 10000-1];	# column element numbers
set level 0
for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
	for {set level 1} {$level <=$NStory} {incr level 1} {
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
			set elemID [expr $N0col  +$level*$Dlevel + $frame*$Dframe+$pier]
			set nodeI [expr  $level*$Dlevel + $frame*$Dframe+$pier]
			set nodeJ  [expr  ($level+1)*$Dlevel + $frame*$Dframe+$pier]
			element nonlinearBeamColumn $elemID $nodeI $nodeJ $numIntgrPts $ColSecTag $IDColTransf;		# columns
		}
	}
}


# beams -- parallel to X-axis
set N0beam 1000000;	# beam element numbers
for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
	for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
		for {set bay 1} {$bay <= $NBay} {incr bay 1} {
			set elemID [expr $N0beam +$level*$Dlevel + $frame*$Dframe+ $bay]
			set nodeI [expr $level*$Dlevel + $frame*$Dframe+ $bay]
			set nodeJ  [expr $level*$Dlevel + $frame*$Dframe+ $bay+1]
			element nonlinearBeamColumn $elemID $nodeI $nodeJ $numIntgrPts $BeamSecTag $IDBeamTransf;	# beams
		}
	}
}

# girders -- parallel to Z-axis
set N0gird 2000000;	# gird element numbers
for {set frame 1} {$frame <=[expr $NFrame-1]} {incr frame 1} {
	for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
		for {set bay 1} {$bay <= $NBay+1} {incr bay 1} {
			set elemID [expr $N0gird + $level*$Dlevel +$frame*$Dframe+ $bay]
			set nodeI [expr   $level*$Dlevel + $frame*$Dframe+ $bay]
			set nodeJ  [expr  $level*$Dlevel + ($frame+1)*$Dframe+ $bay]
			element nonlinearBeamColumn $elemID $nodeI $nodeJ $numIntgrPts $GirdSecTag $IDGirdTransf;		# Girds
		}
	}
}

# --------------------------------------------------------------------------------------------------------------------------------
# Define GRAVITY LOADS, weight and masses
# calculate dead load of frame, assume this to be an internal frame (do LL in a similar manner)
# calculate distributed weight along the beam length
#set GammaConcrete [expr 150*$pcf];   		# Reinforced-Concrete floor slabs, defined above
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
set iFloorWeight ""
set WeightTotal 0.0
set sumWiHi 0.0;		# sum of storey weight times height, for lateral-load distribution

for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
	if {$frame == 1 || $frame == $NFrame}  {
		set GirdWeightFact 1;		# 1x1/2girder on exterior frames
	} else {
		set GirdWeightFact 2;		# 2x1/2girder on interior frames
	}
	for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} { ;		
		set FloorWeight 0.0
		if {$level == [expr $NStory+1]}  {
			set ColWeightFact 1;		# one column in top story
		} else {
			set ColWeightFact 2;		# two columns elsewhere
		}
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {;
			if {$pier == 1 || $pier == [expr $NBay+1]} {
				set BeamWeightFact 1;	# one beam at exterior nodes
			} else {;
				set BeamWeightFact 2;	# two beams elewhere
			}
			set WeightNode [expr $ColWeightFact*$WeightCol/2 + $BeamWeightFact*$WeightBeam/2 + $GirdWeightFact*$WeightGird/2]
			set MassNode [expr $WeightNode/$g];
			set nodeID [expr $level*$Dlevel+$frame*$Dframe+$pier]
			mass $nodeID $MassNode 0. $MassNode 0. 0. 0.;			# define mass
			set FloorWeight [expr $FloorWeight+$WeightNode];
		}
		lappend iFloorWeight $FloorWeight
		set WeightTotal [expr $WeightTotal+ $FloorWeight]
		set sumWiHi [expr $sumWiHi+$FloorWeight*($level-1)*$LCol];		# sum of storey weight times height, for lateral-load distribution
	}
}
set MassTotal [expr $WeightTotal/$g];						# total mass

# --------------------------------------------------------------------------------------------------------------------------------
# LATERAL-LOAD distribution for static pushover analysis
# calculate distribution of lateral load based on mass/weight distributions along building height
# Fj = WjHj/sum(WiHi)  * Weight   at each floor j
set iFj ""
for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} { ;	
	set FloorWeight [lindex $iFloorWeight [expr $level-1-1]];
	set FloorHeight [expr ($level-1)*$LCol];
	lappend iFj [expr $FloorWeight*$FloorHeight/$sumWiHi*$WeightTotal];		# per floor
}
set iNodePush $iMasterNode;		# nodes for pushover/cyclic, vectorized
set iFPush $iFj;				# lateral load for pushover, vectorized

# Define RECORDERS -------------------------------------------------------------
set FreeNodeID [expr $NFrame*$Dframe+($NStory+1)*$Dlevel+($NBay+1)];					# ID: free node
set SupportNodeFirst [lindex $iSupportNode 0];						# ID: first support node
set SupportNodeLast [lindex $iSupportNode [expr [llength $iSupportNode]-1]];			# ID: last support node
set FirstColumn [expr $N0col  + 1*$Dframe+1*$Dlevel +1];							# ID: first column
recorder Node -file $dataDir/DFree.out -time -node $FreeNodeID  -dof 1 2 3 disp;				# displacements of free node
recorder Node -file $dataDir/DBase.out -time -nodeRange $SupportNodeFirst $SupportNodeLast -dof 1 2 3 disp;	# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -nodeRange $SupportNodeFirst $SupportNodeLast -dof 1 2 3 reaction;	# support reaction
recorder Drift -file $dataDir/DrNode.out -time -iNode $SupportNodeFirst  -jNode $FreeNodeID   -dof 1 -perpDirn 2;	# lateral drift
recorder Element -file $dataDir/Fel1.out -time -ele $FirstColumn localForce;					# element forces in local coordinates
recorder Element -file $dataDir/ForceEle1sec1.out -time -ele $FirstColumn section 1 force;			# section forces, axial and moment, node i
recorder Element -file $dataDir/DefoEle1sec1.out -time -ele $FirstColumn section 1 deformation;			# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceEle1sec$numIntgrPts.out -time -ele $FirstColumn section $numIntgrPts force;			# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoEle1sec$numIntgrPts.out -time -ele $FirstColumn section $numIntgrPts deformation;		# section deformations, axial and curvature, node j
set yFiber [expr $HCol/2-$cover];								# fiber location for stress-strain recorder, local coords
set zFiber [expr $BCol/2-$cover];								# fiber location for stress-strain recorder, local coords
recorder Element -file $dataDir/SSconcEle1sec1.out -time -ele $FirstColumn section $numIntgrPts fiber $yFiber $zFiber $IDconcCore  stressStrain;	# steel fiber stress-strain, node i
recorder Element -file $dataDir/SSreinfEle1sec1.out -time -ele $FirstColumn section $numIntgrPts fiber $yFiber $zFiber $IDSteel  stressStrain;	# steel fiber stress-strain, node i

# Define DISPLAY -------------------------------------------------------------
#DisplayModel3D DeformedShape ;	 # options: DeformedShape NodeNumbers ModeShape

# Modifications

# inelastic section for the infill wall

set EminfM 500;    #masonry modulus of elasticity
set sectioninf 10;

section fiberSec $sectioninf -GJ 1e8 {
	  set infmattag 11;
	  set fyfibinf 0.756;
	  set areafibinf 2.214145;
	  set zfibinf 17.967402;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;
	  
	  set infmattag 12;
	  set fyfibinf 0.611;
	  set areafibinf 5.260807;
	  set zfibinf 9.367025;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;
	  	
	  set infmattag 13;
	  set fyfibinf 0.545;
	  set areafibinf 8.324044;
	  set zfibinf 6.63148;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;

	  set infmattag 14;
	  set fyfibinf 0.49;
	  set areafibinf 12.791111;
	  set zfibinf 4.799369;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;

	  set infmattag 15;
	  set fyfibinf 0.396;
	  set areafibinf 30.176721;
	  set zfibinf 2.515478;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;  
	  	  	 
	  set infmattag 16;
	  set fyfibinf 0.396;
	  set areafibinf 30.176721;
	  set zfibinf -2.515478;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;
  	 
	  set infmattag 17;
	  set fyfibinf 0.49;
	  set areafibinf 12.791111;
	  set zfibinf -4.799369;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;

  	  set infmattag 18;
	  set fyfibinf 0.545;
	  set areafibinf 8.324044;
	  set zfibinf -6.63148;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;  	  

  	  set infmattag 19;
	  set fyfibinf 0.611;
	  set areafibinf 5.260807;
	  set zfibinf -9.367025;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;

  	  set infmattag 20;
	  set fyfibinf 0.756;
	  set areafibinf 2.214145;
	  set zfibinf -17.967402;
	  uniaxialMaterial Steel01 $infmattag $fyfibinf $EminfM 0.02;
	  fiber 0.0 $zfibinf $areafibinf $infmattag;

# this fiber in y directtion with very snall area is needed to supply a very small in plane stiffness	  	  
 	  set infmattag 21;
 	  uniaxialMaterial Steel01 $infmattag 1.e40 $EminfM 0.02;
 	  layer straight $infmattag 1 0.0001 1.0 0.0 1.0 0.0;	  	
}
  
# aggregate with rigid torsion stiffness

set sectioninfT 11;
set Torsionmat 90;
uniaxialMaterial Elastic $Torsionmat [expr 1.e12]
section Aggregator $sectioninfT $Torsionmat T -section $sectioninf

# elastic section for pin connection
set sectionpin 12;
set AreainfM 117.5337;
section Elastic $sectionpin $EminfM $AreainfM 1e-3 1e-3 [expr $EminfM/2.5] 1e-3;

# generate the diagonal elements and the midspan nodes placing the infill walls in the exterior frames in XY plane
puts "ok"
#midspan nodes and the masses at the midspan nodes
set massinf  [expr 5.221*$kip/$g];
set Dlevel 10000;
set Dframe 100;	
foreach frame {1 $NFrame}  {
	set Z [expr ($frame-1)*$LGird];
	for {set level 1} {$level <=[expr $NStory]} {incr level 1} {
		set Y [expr ($level-1+0.5)*$LCol];
		for {set pier 1} {$pier <= [expr $NBay]} {incr pier 1} {
			set X [expr ($pier-1+0.5)*$LBeam];
			set nodeID [expr $level*$Dlevel+1000+$frame*$Dframe+$pier]; 
			node $nodeID $X $Y $Z -mass 0.0 0.0 [expr $massinf] 0.0 0.0 0.0; 	
		}
	}
}

# diagonal elements and the recorders for removal
set InertiainfM 4056.031;
set rinfM [expr pow((pow($LCol,2.0)+pow($LBeam,2.0)),0.5)];
set InfillTransf 9;
geomTransf Linear $InfillTransf 0 0 -1
set infnum 0;
set infnum2 200;
set fileremoval "Dispwall1-cg.tcl"; 
foreach frame {1 $NFrame}  {
	for {set level 1} {$level <=$NStory} {incr level 1} {
		for {set pier 1} {$pier <= [expr $NBay]} {incr pier 1} {             
			set nodeI [expr  $level*$Dlevel + $frame*$Dframe+$pier]
			set nodeCM [expr  $level*$Dlevel + 1000+ $frame*$Dframe+$pier]
			set nodeJ  [expr  ($level+1)*$Dlevel + $frame*$Dframe+$pier+1]
			set infnum [expr $infnum+1]
            element beamWithHinges $infnum $nodeCM  $nodeI $sectioninfT [expr $rinfM*0.1] $sectionpin [expr $rinfM*0.05] $EminfM $AreainfM 1.e-5 $InertiainfM   [expr $EminfM/2.5] $Ubig $InfillTransf
            set infnum2 [expr $infnum2+1]
            element beamWithHinges $infnum2 $nodeCM  $nodeJ $sectioninfT [expr $rinfM*0.1] $sectionpin [expr $rinfM*0.05] $EminfM $AreainfM 1.e-5 $InertiainfM   [expr $EminfM/2.5] $Ubig $InfillTransf

#         recorders for collapse removal			
      	    recorder Collapse -ele $infnum   -time  -crit INFILLWALL  -file $dataDir/CollapseSequence.out  -file_infill $fileremoval -global_gravaxis 2 -checknodes $nodeI $nodeCM $nodeJ
      	    recorder Collapse -ele $infnum2  -time  -crit INFILLWALL   -file_infill $fileremoval -global_gravaxis 2 -checknodes $nodeI $nodeCM $nodeJ 
			recorder Collapse -ele $infnum $infnum2 -node $nodeCM
#          recorders for the displacements of the midspan nodes 
            recorder Node -file $dataDir/Midspan-$nodeCM.out -time -node $nodeCM -dof 1 2 3 disp;			
		}
	}
}

# end of modifications

# GRAVITY -------------------------------------------------------------
# define GRAVITY load applied to beams and columns -- eleLoad applies loads in local coordinate axis
pattern Plain 101 Linear {
	for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
		for {set level 1} {$level <=$NStory} {incr level 1} {
			for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
				set elemID [expr $N0col  + $level*$Dlevel +$frame*$Dframe+$pier]
					eleLoad -ele $elemID -type -beamUniform 0. 0. -$QdlCol; 	# COLUMNS		}
		}
	}
	for {set frame 1} {$frame <=[expr $NFrame]} {incr frame 1} {
		for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
			for {set bay 1} {$bay <= $NBay} {incr bay 1} {
				set elemID [expr $N0beam + $level*$Dlevel +$frame*$Dframe+ $bay]
					eleLoad -ele $elemID  -type -beamUniform -$QdlBeam 0.; 	# BEAMS
			}
		}
	}
	for {set frame 1} {$frame <=[expr $NFrame-1]} {incr frame 1} {
		for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
			for {set bay 1} {$bay <= $NBay+1} {incr bay 1} {
				set elemID [expr $N0gird + $level*$Dlevel +$frame*$Dframe+ $bay]
				eleLoad -ele $elemID  -type -beamUniform -$QdlGird 0.;	# GIRDS
			}
		}
	}

}
puts goGravity
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

# -------------------------------------------------------------
puts "Model Built"

#
