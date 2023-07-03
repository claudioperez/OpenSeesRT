# --------------------------------------------------------------------------------------------------
# Example: 2-Story Steel Moment Frame with Concentrated Plasticity
# Centerline Model with Concentrated Plastic Hinges at Beam-Column Joint
# Created by:  Laura Eads, Stanford University, 2010
# Units: kips, inches, seconds

# Element and Node ID conventions:
#	1xy = frame columns with springs at both ends
#	2xy = frame beams with springs at both ends
#	6xy = trusses linking frame and P-delta column
#	7xy = P-delta columns
#	3,xya = frame column rotational springs
#	4,xya = frame beam rotational springs
#	5,xya = P-delta column rotational springs
#	where:
#		x = Pier or Bay #
#		y = Floor or Story #
#		a = an integer describing the location relative to beam-column joint (see description where elements and nodes are defined)

###################################################################################################
#          Set Up & Source Definition									  
###################################################################################################
	wipe all;							# clear memory of past model definitions
	model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm = #dimension, ndf = #dofs
	source DisplayModel2D.tcl;			# procedure for displaying a 2D perspective of model
	source DisplayPlane.tcl;			# procedure for displaying a plane in a model
	source rotSpring2DModIKModel.tcl;	# procedure for defining a rotational spring (zero-length element)
	source rotLeaningCol.tcl;			# procedure for defining a rotational spring (zero-length element) with very small stiffness
	
###################################################################################################
#          Define Analysis Type										  
###################################################################################################
# Define type of analysis:  "pushover" = pushover
	set analysisType "pushover";
	if {$analysisType == "pushover"} {
		set dataDir Concentrated-Pushover-Output;	# name of output folder
		file mkdir $dataDir;						# create output folder
	}
	
###################################################################################################
#          Define Building Geometry, Nodes, and Constraints											  
###################################################################################################
# define structure-geometry parameters
	set NStories 2;						# number of stories
	set NBays 1;						# number of frame bays (excludes bay for P-delta column)
	set WBay      [expr 30.0*12.0];		# bay width in inches
	set HStory1   [expr 15.0*12.0];		# 1st story height in inches
	set HStoryTyp [expr 12.0*12.0];		# story height of other stories in inches
	set HBuilding [expr $HStory1 + ($NStories-1)*$HStoryTyp];	# height of building

# calculate locations of beam/column joints:
	set Pier1  0.0;		# leftmost column line
	set Pier2  [expr $Pier1 + $WBay];
	set Pier3  [expr $Pier2 + $WBay];	# P-delta column line	
	set Floor1 0.0;		# ground floor
	set Floor2 [expr $Floor1 + $HStory1];
	set Floor3 [expr $Floor2 + $HStoryTyp];

# calculate joint offset distance for beam plastic hinges
	set phlat23 [expr 0.0];		# lateral dist from beam-col joint to loc of hinge on Floor 2

# calculate nodal masses -- lump floor masses at frame nodes
	set g 386.2;				# acceleration due to gravity
	set Floor2Weight 535.0;		# weight of Floor 2 in kips
	set Floor3Weight 525.0;		# weight of Floor 3 in kips
	set WBuilding  [expr $Floor2Weight + $Floor3Weight];# total building weight
	set NodalMass2 [expr ($Floor2Weight/$g) / (2.0)];	# mass at each node on Floor 2
	set NodalMass3 [expr ($Floor3Weight/$g) / (2.0)];	# mass at each node on Floor 3
	set Negligible 1e-9;	# a very smnumber to avoid problems with zero

# define nodes and assign masses to beam-column intersections of frame
	# command:  node nodeID xcoord ycoord -mass mass_dof1 mass_dof2 mass_dof3
	# nodeID convention:  "xy" where x = Pier # and y = Floor # 
	node 11 $Pier1 $Floor1;
	node 21 $Pier2 $Floor1;
	node 31 $Pier3 $Floor1;
	node 12 $Pier1 $Floor2 -mass $NodalMass2 $Negligible $Negligible;
	node 22 $Pier2 $Floor2 -mass $NodalMass2 $Negligible $Negligible;
	node 32 $Pier3 $Floor2;
	node 13 $Pier1 $Floor3 -mass $NodalMass3 $Negligible $Negligible;
	node 23 $Pier2 $Floor3 -mass $NodalMass3 $Negligible $Negligible;
	node 33 $Pier3 $Floor3;

# define extra nodes for plastic hinge rotational springs
	# nodeID convention:  "xya" where x = Pier #, y = Floor #, a = location relative to beam-column joint
	# "a" convention: 2 = left; 3 = right;
	# "a" convention: 6 = below; 7 = above; 
	# column hinges at bottom of Story 1 (base)
	node 117 $Pier1 $Floor1;
	node 217 $Pier2 $Floor1;
	# column hinges at top of Story 1
	node 126 $Pier1 $Floor2;
	node 226 $Pier2 $Floor2;
	node 326 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	# column hinges at bottom of Story 2
	node 127 $Pier1 $Floor2;
	node 227 $Pier2 $Floor2;
	node 327 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	# column hinges at top of Story 2
	node 136 $Pier1 $Floor3;
	node 236 $Pier2 $Floor3;
	node 336 $Pier3 $Floor3;	# zero-stiffness spring will be used on p-delta column

	# beam hinges at Floor 2
	node 122 [expr $Pier1 + $phlat23] $Floor2;
	node 223 [expr $Pier2 - $phlat23] $Floor2;
	# beam hinges at Floor 3
	node 132 [expr $Pier1 + $phlat23] $Floor3;
	node 233 [expr $Pier2 - $phlat23] $Floor3;

# constrain beam-column joints in a floor to have the same lateral displacement using the "equalDOF" command
	# command: equalDOF $MasterNodeID $SlaveNodeID $dof1 $dof2...
	set dof1 1;	# constrain movement in dof 1 (x-direction)
	equalDOF 12 22 $dof1;	# Floor 2:  Pier 1 to Pier 2
	equalDOF 12 32 $dof1;	# Floor 2:  Pier 1 to Pier 3
	equalDOF 13 23 $dof1;	# Floor 3:  Pier 1 to Pier 2
	equalDOF 13 33 $dof1;	# Floor 3:  Pier 1 to Pier 3

# assign boundary condidtions 
	# command:  fix nodeID dxFixity dyFixity rzFixity
	# fixity values: 1 = constrained; 0 = unconstrained
	# fix the base of the building; pin P-delta column at base
	fix 11 1 1 1;
	fix 21 1 1 1;
	fix 31 1 1 0;	# P-delta column is pinned

###################################################################################################
#          Define Section Properties and Elements													  
###################################################################################################
# define material properties
	set Es 29000.0;			# steel Young's modulus

# define column section W24x131 for Story 1 & 2
	set Acol_12  38.5;		# cross-sectional area
	set Icol_12  4020.0;	# moment of inertia
	set Mycol_12 20350.0;	# yield moment

# define beam section W27x102 for Floor 2 & 3
	set Abeam_23  30.0;		# cross-sectional area (full section properties)
	set Ibeam_23  3620.0;	# moment of inertia  (full section properties)
	set Mybeam_23 10938.0;	# yield moment at plastic hinge location (i.e., My of RBS section, if used)
	# note: In this example the hinges form right at the beam-column joint, so using an RBS doesn't make sense; 
	#		however, it is done here simply for illustrative purposes.
	
# determine stiffness modifications to equate the stiffness of the spring-elastic element-spring subassembly to the stiffness of the actual frame member
	# Reference:  Ibarra, L. F., and Krawinkler, H. (2005). "Global collapse of frame structures under seismic excitations," Technical Report 152,
	#             The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
	# calculate modified section properties to account for spring stiffness being in series with the elastic element stiffness
	set n 10.0;		# stiffness multiplier for rotational spring

	# calculate modified moment of inertia for elastic elements
	set Icol_12mod  [expr $Icol_12*($n+1.0)/$n];	# modified moment of inertia for columns in Story 1 & 2
	set Ibeam_23mod [expr $Ibeam_23*($n+1.0)/$n];	# modified moment of inertia for beams in Floor 2 & 3
	# calculate modified rotational stiffness for plastic hinge springs
	set Ks_col_1   [expr $n*6.0*$Es*$Icol_12mod/$HStory1];		# rotational stiffness of Story 1 column springs
	set Ks_col_2   [expr $n*6.0*$Es*$Icol_12mod/$HStoryTyp];	# rotational stiffness of Story 2 column springs
	set Ks_beam_23 [expr $n*6.0*$Es*$Ibeam_23mod/$WBay];		# rotational stiffness of Floor 2 & 3 beam springs
	
# set up geometric transformations of element
	set PDeltaTransf 1;
	geomTransf PDelta $PDeltaTransf; 	# PDelta transformation

# define elastic column elements using "element" command
	# command: element elasticBeamColumn $eleID $iNode $jNode $A $E $I $transfID
	# eleID convention:  "1xy" where 1 = col, x = Pier #, y = Story #
	# Columns Story 1
	element elasticBeamColumn  111  117 126 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 1
	element elasticBeamColumn  121  217 226 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 2
	# Columns Story 2
	element elasticBeamColumn  112  127 136 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 1
	element elasticBeamColumn  122  227 236 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 2
	
# define elastic beam elements
	# eleID convention:  "2xy" where 2 = beam, x = Bay #, y = Floor #
	# Beams Story 1
	element elasticBeamColumn  212  122 223 $Abeam_23 $Es $Ibeam_23mod $PDeltaTransf;
	# Beams Story 2
	element elasticBeamColumn  222  132 233 $Abeam_23 $Es $Ibeam_23mod $PDeltaTransf;
	
# define p-delta columns and rigid links
	set TrussMatID 600;		# define a material ID
	set Arigid 1000.0;		# define area of truss section (make much larger than A of frame elements)
	set Irigid 100000.0;	# moment of inertia for p-delta columns  (make much larger than I of frame elements)
	uniaxialMaterial Elastic $TrussMatID $Es;		# define truss material
	# rigid links
	# command: element truss $eleID $iNode $jNode $A $materialID
	# eleID convention:  6xy, 6 = truss link, x = Bay #, y = Floor #
	element truss  622 22 32 $Arigid $TrussMatID;	# Floor 2
	element truss  623 23 33 $Arigid $TrussMatID;	# Floor 3
	
	# p-delta columns
	# eleID convention:  7xy, 7 = p-delta columns, x = Pier #, y = Story #
	element elasticBeamColumn  731  31  326 $Arigid $Es $Irigid $PDeltaTransf;	# Story 1
	element elasticBeamColumn  732  327 336 $Arigid $Es $Irigid $PDeltaTransf;	# Story 2
	
# display the model with the node numbers
	DisplayModel2D NodeNumbers
	
###################################################################################################
#          Define Rotational Springs for Plastic Hinges												  
###################################################################################################
# define rotational spring properties and create spring elements using "rotSpring2DModIKModel" procedure
	# rotSpring2DModIKModel creates a uniaxial material spring with a bilinear response based on Modified Ibarra Krawinkler Deterioration Model
	# references provided in rotSpring2DModIKModel.tcl
	# input values for Story 1 column springs
	set McMy 1.05;			# ratio of capping moment to yield moment, Mc / My
	set LS 1000.0;			# basic strength deterioration (a very large # = no cyclic deterioration)
	set LK 1000.0;			# unloading stiffness deterioration (a very large # = no cyclic deterioration)
	set LA 1000.0;			# accelerated reloading stiffness deterioration (a very large # = no cyclic deterioration)
	set LD 1000.0;			# post-capping strength deterioration (a very large # = no deterioration)
	set cS 1.0;				# exponent for basic strength deterioration (c = 1.0 for no deterioration)
	set cK 1.0;				# exponent for unloading stiffness deterioration (c = 1.0 for no deterioration)
	set cA 1.0;				# exponent for accelerated reloading stiffness deterioration (c = 1.0 for no deterioration)
	set cD 1.0;				# exponent for post-capping strength deterioration (c = 1.0 for no deterioration)
	set th_pP 0.025;		# plastic rot capacity for pos loading
	set th_pN 0.025;		# plastic rot capacity for neg loading
	set th_pcP 0.3;			# post-capping rot capacity for pos loading
	set th_pcN 0.3;			# post-capping rot capacity for neg loading
	set ResP 0.4;			# residual strength ratio for pos loading
	set ResN 0.4;			# residual strength ratio for neg loading
	set th_uP 0.4;			# ultimate rot capacity for pos loading
	set th_uN 0.4;			# ultimate rot capacity for neg loading
	set DP 1.0;				# rate of cyclic deterioration for pos loading
	set DN 1.0;				# rate of cyclic deterioration for neg loading
	set a_mem [expr ($n+1.0)*($Mycol_12*($McMy-1.0)) / ($Ks_col_1*$th_pP)];	# strain hardening ratio of spring
	set b [expr ($a_mem)/(1.0+$n*(1.0-$a_mem))];							# modified strain hardening ratio of spring (Ibarra & Krawinkler 2005, note: Eqn B.5 is incorrect)

# define column springs
	# Spring ID: "3xya", where 3 = col spring, x = Pier #, y = Story #, a = location in story
	# "a" convention: 1 = bottom of story, 2 = top of story
	# command: rotSpring2DModIKModel	id    ndR  ndC     K   asPos  asNeg  MyPos      MyNeg      LS    LK    LA    LD   cS   cK   cA   cD  th_p+   th_p-   th_pc+   th_pc-  Res+   Res-   th_u+   th_u-    D+     D-
	# col springs @ bottom of Story 1 (at base)
	rotSpring2DModIKModel 3111 11 117 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 3211 21 217 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#col springs @ top of Story 1 (below Floor 2)
	rotSpring2DModIKModel 3112 12 126 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 3212 22 226 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

	# recompute strain hardening since Story 2 is not the same height as Story 1
	set a_mem [expr ($n+1.0)*($Mycol_12*($McMy-1.0)) / ($Ks_col_2*$th_pP)];	# strain hardening ratio of spring
	set b [expr ($a_mem)/(1.0+$n*(1.0-$a_mem))];							# modified strain hardening ratio of spring (Ibarra & Krawinkler 2005, note: there is mistake in Eqn B.5)
	# col springs @ bottom of Story 2 (above Floor 2)
	rotSpring2DModIKModel 3121 12 127 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 3221 22 227 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#col springs @ top of Story 2 (below Floor 3)
	rotSpring2DModIKModel 3122 13 136 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 3222 23 236 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	
	# create region for frame column springs
	# command: region $regionID -ele $ele_1_ID $ele_2_ID...
	region 1 -ele 3111 3211 3112 3212 3121 3221 3122 3222;
	
# define beam springs
	# Spring ID: "4xya", where 4 = beam spring, x = Bay #, y = Floor #, a = location in bay
	# "a" convention: 1 = left end, 2 = right end
	# redefine the rotations since they are not the same
	set th_pP 0.02;
	set th_pN 0.02;
	set th_pcP 0.16;
	set th_pcN 0.16;
	set a_mem [expr ($n+1.0)*($Mybeam_23*($McMy-1.0)) / ($Ks_beam_23*$th_pP)];	# strain hardening ratio of spring
	set b [expr ($a_mem)/(1.0+$n*(1.0-$a_mem))];								# modified strain hardening ratio of spring (Ibarra & Krawinkler 2005, note: there is mistake in Eqn B.5)
	#beam springs at Floor 2
	rotSpring2DModIKModel 4121 12 122 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 4122 22 223 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#beam springs at Floor 3
	rotSpring2DModIKModel 4131 13 132 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	rotSpring2DModIKModel 4132 23 233 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	
	# create region for beam springs
	region 2 -ele 4121 4122 4131 4132;
	
# define p-delta column spring: zero-stiffness elastic spring	
	#Spring ID: "5xya" where 5 = leaning column spring, x = Pier #, y = Story #, a = location in story
	# "a" convention: 1 = bottom of story, 2 = top of story
	# rotLeaningCol ElemID ndR ndC 
	rotLeaningCol 5312 32 326;	# top of Story 1
	rotLeaningCol 5321 32 327;	# bottom of Story 2
	rotLeaningCol 5322 33 336;	# top of Story 2
	
	# create region for P-Delta column springs
	region 3 -ele 5312 5321 5322;
	
############################################################################
#                       Eigenvalue Analysis                    			   
############################################################################
	set pi [expr 2.0*asin(1.0)];						# Definition of pi
	set nEigenI 1;										# mode i = 1
	set nEigenJ 2;										# mode j = 2
	set lambdaN [eigen [expr $nEigenJ]];				# eigenvalue analysis for nEigenJ modes
	set lambdaI [lindex $lambdaN [expr 0]];				# eigenvalue mode i = 1
	set lambdaJ [lindex $lambdaN [expr $nEigenJ-1]];	# eigenvalue mode j = 2
	set w1 [expr pow($lambdaI,0.5)];					# w1 (1st mode circular frequency)
	set w2 [expr pow($lambdaJ,0.5)];					# w2 (2nd mode circular frequency)
	set T1 [expr 2.0*$pi/$w1];							# 1st mode period of the structure
	set T2 [expr 2.0*$pi/$w2];							# 2nd mode period of the structure
	puts "T1 = $T1 s";									# display the first mode period in the command window
	puts "T2 = $T2 s";									# display the second mode period in the command window
	
############################################################################
#              Gravity Loads & Gravity Analysis
############################################################################
# apply gravity loads
	#command: pattern PatternType $PatternID TimeSeriesType
	pattern Plain 101 Constant {
		
		# point loads on leaning column nodes
		# command: load node Fx Fy Mz
		set P_PD2 [expr -398.02];	# Floor 2
		set P_PD3 [expr -391.31];	# Floor 3
		load 32 0.0 $P_PD2 0.0;		# Floor 2
		load 33 0.0 $P_PD3 0.0;		# Floor 3
		
		# point loads on frame column nodes
		set P_F2 [expr 0.5*(-1.0*$Floor2Weight-$P_PD2)];	# load on each frame node in Floor 2
		set P_F3 [expr 0.5*(-1.0*$Floor3Weight-$P_PD3)];	# load on each frame node in Floor 3
		# Floor 2 loads
		load 12 0.0 $P_F2 0.0;
		load 22 0.0 $P_F2 0.0;		
		# Floor 3 loads		
		load 13 0.0 $P_F3 0.0;
		load 23 0.0 $P_F3 0.0;
	}

# Gravity-analysis: load-controlled static analysis
	set Tol 1.0e-6;							# convergence tolerance for test
	constraints Plain;						# how it handles boundary conditions
	numberer RCM;							# renumber dof's to minimize band-width (optimization)
	system BandGeneral;						# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	test NormDispIncr $Tol 6;				# determine if convergence has been achieved at the end of an iteration step
	algorithm Newton;						# use Newton's solution algorithm: updates tangent stiffness at every iteration
	set NstepGravity 10;					# apply gravity in 10 steps
	set DGravity [expr 1.0/$NstepGravity];	# load increment
	integrator LoadControl $DGravity;		# determine the next time step for an analysis
	analysis Static;						# define type of analysis static or transient
	analyze $NstepGravity;					# apply gravity

	# maintain constant gravity loads and reset time to zero
	loadConst -time 0.0
	puts "Model Built"
	
############################################################################
#              Recorders					                			   
############################################################################
# record drift histories
	# drift recorder command: recorder Drift -file $filename -iNode $NodeI_ID -jNode $NodeJ_ID -dof $dof -perpDirn $Record.drift.perpendicular.to.this.direction
	recorder Drift -file $dataDir/Drift-Story1.out -iNode 11 -jNode 12 -dof 1 -perpDirn 2;
	recorder Drift -file $dataDir/Drift-Story2.out -iNode 12 -jNode 13 -dof 1 -perpDirn 2;
	recorder Drift -file $dataDir/Drift-Roof.out -iNode 11 -jNode 13 -dof 1 -perpDirn 2;
	
# record base shear reactions
	recorder Node -file $dataDir/Vbase.out -node 117 217 31 -dof 1 reaction;
	
# record story 1 column forces in global coordinates 
	recorder Element -file $dataDir/Fcol111.out -ele 111 force;
	recorder Element -file $dataDir/Fcol121.out -ele 121 force;
	recorder Element -file $dataDir/Fcol731.out -ele 731 force;
	
# record response history of all frame column springs (one file for moment, one for rotation)
	recorder Element -file $dataDir/MRFcol-Mom-Hist.out -region 1 force;
	recorder Element -file $dataDir/MRFcol-Rot-Hist.out -region 1 deformation;
	
# record response history of all frame beam springs (one file for moment, one for rotation)
	recorder Element -file $dataDir/MRFbeam-Mom-Hist.out -region 2 force;
	recorder Element -file $dataDir/MRFbeam-Rot-Hist.out -region 2 deformation;
	
#######################################################################################
#                                                                                     #
#                              Analysis Section			                              #
#                                                                                     #
#######################################################################################

############################################################################
#              Pushover Analysis                			   			   #
############################################################################
if {$analysisType == "pushover"} { 
	puts "Running Pushover..."
# assign lateral loads and create load pattern:  use ASCE 7-10 distribution
	set lat2 16.255;	# force on each frame node in Floor 2
	set lat3 31.636;	# force on each frame node in Floor 3
	pattern Plain 200 Linear {			
					load 12 $lat2 0.0 0.0;
					load 22 $lat2 0.0 0.0;
					load 13 $lat3 0.0 0.0;
					load 23 $lat3 0.0 0.0;
	}
	
# display deformed shape:
	set ViewScale 5;
	DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model

# displacement parameters
	set IDctrlNode 13;					# node where disp is read for disp control
	set IDctrlDOF 1;					# degree of freedom read for disp control (1 = x displacement)
	set Dmax [expr 0.1*$HBuilding];		# maximum displacement of pushover: 10% roof drift
	set Dincr [expr 0.01];				# displacement increment

# analysis commands
	constraints Plain;					# how it handles boundary conditions
	numberer RCM;						# renumber dof's to minimize band-width (optimization)
	system BandGeneral;					# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	test NormUnbalance 1.0e-6 400;		# tolerance, max iterations
	algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
	integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr;	# use displacement-controlled analysis
	analysis Static;					# define type of analysis: static for pushover
	set Nsteps [expr int($Dmax/$Dincr)];# number of pushover analysis steps
	set ok [analyze $Nsteps];			# this will return zero if no convergence problems were encountered
	puts "Pushover complete";			# display this message in the command window
} 	
	
	
wipe all;