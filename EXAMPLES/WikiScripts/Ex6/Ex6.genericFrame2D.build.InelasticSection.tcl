# --------------------------------------------------------------------------------------------------
# Example 6. 2D Frame --  Inelastic Uniaxial Section
#		Silvia Mazzoni & Frank McKenna, 2006
# nonlinearBeamColumn element, inelastic fiber section
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Output;			# set up name of data directory (can remove this)
file mkdir $dataDir; 			# create data directory
set GMdir "Motions";			# ground-motion file directory
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source LibUnits.tcl;			# define units

# define GEOMETRY -------------------------------------------------------------
# define structure-geometry paramters
set LCol [expr 14*$ft];		# column height
set LBeam [expr 24*$ft];		# beam length
set NStory 3;			# number of stories above ground level -------------- you can change this.
set NBay 3;			# number of bays (max 9) ------------------------------you can change this.

# define NODAL COORDINATES
for {set level 1} {$level <=[expr $NStory+1]} {incr level 1} {
	set Y [expr ($level-1)*$LCol];
	for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
		set X [expr ($pier-1)*$LBeam];
		set nodeID [expr $level*10+$pier]
		node $nodeID $X $Y;		# actually define node
	}
}
# determine support nodes where ground motions are input, for multiple-support excitation
set iSupportNode ""
set level 1
for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
	set nodeID [expr $level*10+$pier]
	lappend iSupportNode $nodeID
}



# BOUNDARY CONDITIONS
fixY 0.0 1 1 0;		# pin all Y=0.0 nodes

# calculated MODEL PARAMETERS, particular to this model
puts "Number of Stories: $NStory Number of bays: $NBay"
# Set up parameters that are particular to the model for displacement control
set IDctrlNode [expr ($NStory+1)*10+1];		# node where displacement is read for displacement control
set IDctrlDOF 1;		# degree of freedom of displacement read for displacement control
set LBuilding [expr $NStory*$LCol];	# total building height

# Define ELEMENTS & SECTIONS  -------------------------------------------------------------
set ColSecTag 1;				# assign a tag number to the column section tag
set ColMatTagFlex 2;			# assign a tag number to the column flexural behavior
set ColMatTagAxial 3;			# assign a tag number to the column axial behavior	
set BeamSecTag 4;				# assign a tag number to the beam section tag
set BeamMatTagFlex 5;			# assign a tag number to the beam flexural behavior
set BeamMatTagAxial 6;			# assign a tag number to the beam axial behavior	

# define MATERIAL properties ----------------------------------------
set Fy [expr 6.0*$ksi]
set Es [expr 29000*$ksi];		# Steel Young's Modulus
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  	# Torsional stiffness Modulus

# COLUMN section W27x114
set AgCol [expr 33.5*pow($in,2)];		# cross-sectional area
set IzCol [expr 4090.*pow($in,4)];		# moment of Inertia
set EICol [expr $Es*$IzCol];				# EI, for moment-curvature relationship
set EACol [expr $Es*$AgCol];				# EA, for axial-force-strain relationship
set MyCol [expr 2e4*$kip*$in];			# yield moment
set PhiYCol [expr 0.25e-3/$in];			# yield curvature
set PhiYCol [expr $MyCol/$EICol];			# yield curvature
set EIColCrack [expr $MyCol/$PhiYCol];		# cracked section inertia
set b 0.01 ;					# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $ColMatTagFlex $MyCol $EIColCrack $b; 		# bilinear behavior for flexure
uniaxialMaterial Elastic $ColMatTagAxial $EACol;				# this is not used as a material, this is an axial-force-strain response
section Aggregator $ColSecTag $ColMatTagAxial P $ColMatTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)

# BEAM section W24x94
set AgBeam [expr 27.7*pow($in,2)];		# cross-sectional area
set IzBeam [expr 2700.*pow($in,4)];		# moment of Inertia
set EIBeam [expr $Es*$IzBeam];			# EI, for moment-curvature relationship
set EABeam [expr $Es*$AgBeam];			# EA, for axial-force-strain relationship
set MyBeam [expr 1.5e4*$kip*$in];			# yield moment
set PhiYBeam [expr 0.25e-3/$in];			# yield curvature
set PhiYBeam [expr $MyBeam/$EIBeam];		# yield curvature
set EIBeamCrack [expr $MyBeam/$PhiYBeam];		# cracked section inertia
set b 0.01 ;					# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $BeamMatTagFlex $MyBeam $EIBeamCrack $b; 		# bilinear behavior for flexure
uniaxialMaterial Elastic $BeamMatTagAxial $EABeam;				# this is not used as a material, this is an axial-force-strain response
section Aggregator $BeamSecTag $BeamMatTagAxial P $BeamMatTagFlex Mz;	# combine axial and flexural behavior into one section (no P-M interaction here)


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
set N0col 100;	# column element numbers
set level 0
for {set level 1} {$level <=$NStory} {incr level 1} {
	for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
		set elemID [expr $N0col  + $level*10 +$pier]
		set nodeI [expr  $level*10 + $pier]
		set nodeJ  [expr  ($level+1)*10 + $pier]
		element nonlinearBeamColumn $elemID $nodeI $nodeJ $np $ColSecTag $IDColTransf;		# columns
	}
}
# beams
set N0beam 200;	# beam element numbers
set M0 0
for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
	for {set bay 1} {$bay <= $NBay} {incr bay 1} {
		set elemID [expr $N0beam + $level*10 +$bay]
		set nodeI [expr  $M0 + $level*10 + $bay]
		set nodeJ  [expr  $M0 + $level*10 + $bay+1]
		element nonlinearBeamColumn $elemID $nodeI $nodeJ $np $BeamSecTag $IDBeamTransf;	# beams
	}
}
                  
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
set iFloorWeight ""
set WeightTotal 0.0
set sumWiHi 0.0;		# sum of storey weight times height, for lateral-load distribution
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
		set WeightNode [expr $ColWeightFact*$WeightCol/2 + $BeamWeightFact*$WeightBeam/2]
		set MassNode [expr $WeightNode/$g];
		set nodeID [expr $level*10+$pier]
		mass $nodeID $MassNode 0.0 0.0; # 0.0 0.0 0.0;			# define mass
		set FloorWeight [expr $FloorWeight+$WeightNode];
	}
	lappend iFloorWeight $FloorWeight
	set WeightTotal [expr $WeightTotal+ $FloorWeight]
	set sumWiHi [expr $sumWiHi+$FloorWeight*($level-1)*$LCol];		# sum of storey weight times height, for lateral-load distribution
}
set MassTotal [expr $WeightTotal/$g];						# total mass

# LATERAL-LOAD distribution for static pushover analysis
# calculate distribution of lateral load based on mass/weight distributions along building height
# Fj = WjHj/sum(WiHi)  * Weight   at each floor j
set iFj ""
for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} { ;	
	set FloorWeight [lindex $iFloorWeight [expr $level-1-1]];
	set FloorHeight [expr ($level-1)*$LCol];
	set NodeFactor [expr $NBay+1];
	lappend iFj [expr $FloorWeight*$FloorHeight/$sumWiHi/$NodeFactor*$WeightTotal];		# per node per floor
}
# create node and load vectors for lateral-load distribution in static analysis
set iFPush ""
set iNodePush ""
for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
	set FPush [lindex $iFj [expr $level-1-1]];		# lateral load coefficient
	for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
		set nodeID [expr $level*10+$pier]
		lappend iNodePush $nodeID
		lappend iFPush $FPush
	}
}

# Define RECORDERS -------------------------------------------------------------
set FreeNodeID [expr ($NStory+1)*10+($NBay+1)];					# ID: free node
set SupportNodeFirst [lindex $iSupportNode 0];						# ID: first support node
set SupportNodeLast [lindex $iSupportNode [expr [llength $iSupportNode]-1]];			# ID: last support node
set FirstColumn [expr $N0col+1*10+1];							# ID: first column
recorder Node -file $dataDir/DFree.out -time -node $FreeNodeID  -dof 1 2 3 disp;				# displacements of free node
recorder Node -file $dataDir/DBase.out -time -nodeRange $SupportNodeFirst $SupportNodeLast -dof 1 2 3 disp;	# displacements of support nodes
recorder Node -file $dataDir/RBase.out -time -nodeRange $SupportNodeFirst $SupportNodeLast -dof 1 2 3 reaction;	# support reaction
recorder Drift -file $dataDir/DrNode.out -time -iNode $SupportNodeFirst  -jNode $FreeNodeID  -dof 1 -perpDirn 2;	# lateral drift
recorder Element -file $dataDir/Fel1.out -time -ele $FirstColumn localForce;					# element forces in local coordinates
recorder Element -file $dataDir/ForceEle1sec1.out -time -ele $FirstColumn section 1 force;			# section forces, axial and moment, node i
recorder Element -file $dataDir/DefoEle1sec1.out -time -ele $FirstColumn section 1 deformation;			# section deformations, axial and curvature, node i
recorder Element -file $dataDir/ForceEle1sec$np.out -time -ele $FirstColumn section $np force;			# section forces, axial and moment, node j
recorder Element -file $dataDir/DefoEle1sec$np.out -time -ele $FirstColumn section $np deformation;		# section deformations, axial and curvature, node j

# Define DISPLAY -------------------------------------------------------------
#DisplayModel2D NodeNumbers

# define GRAVITY -------------------------------------------------------------
# GRAVITY LOADS # define gravity load applied to beams and columns -- eleLoad applies loads in local coordinate axis
pattern Plain 101 Linear {
	for {set level 1} {$level <=$NStory} {incr level 1} {
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
			set elemID [expr $N0col  + $level*10 +$pier]
			eleLoad -ele $elemID -type -beamUniform 0 -$QdlCol; 	# COLUMNS
		}
	}
	for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
		for {set bay 1} {$bay <= $NBay} {incr bay 1} {
			set elemID [expr $N0beam + $level*10 +$bay]
			eleLoad -ele $elemID  -type -beamUniform -$QdlBeam; 	# BEAMS
		}
	}
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


