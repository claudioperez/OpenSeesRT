# --------------------------------------------------------------------------------------------------
# Example: 2-Story Moment Frame with Distributed Plasticity
# Centerline Model with Distributed Plastic Hinges starting at Beam-Column Joint
# Created by:  Laura Eads, Stanford University, 2010
# Units: kips, inches, seconds

# Element and Node ID conventions:
#	1xy = frame columns
#	2xy = frame beams
#   6xy = trusses linking frame and P-delta column
#   7xy = P-delta columns
#	5,xya = P-delta column rotational springs
#	where:
#		x = Pier or Bay #
#   	y = Floor or Story #
#		a = an integer describing the location relative to beam-column joint (see description where elements and nodes are defined)

###################################################################################################
#          Set Up & Source Definition									  
###################################################################################################
	wipe all;							# clear memory of past model definitions
	model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm = #dimension, ndf = #dofs
	source DisplayModel2D.tcl;			# procedure for displaying a 2D perspective of model
	source DisplayPlane.tcl;			# procedure for displaying a plane in the model
	source rotSect2DModIKModel.tcl;		# procedure for defining bilinear plastic hinge section
	source rotLeaningCol.tcl;			# procedure for defining a rotational spring (zero-length element) with very small stiffness
	
###################################################################################################
#          Define Analysis Type										  
###################################################################################################
# Define type of analysis:  "pushover" = pushover
	set analysisType "pushover";
	if {$analysisType == "pushover"} {
		set dataDir Distributed-Pushover-Output;	# name of output folder
		file mkdir $dataDir; 						# create output folder
	}
	
###################################################################################################
#          Define Building Geometry, Nodes, and Constraints											  
###################################################################################################
# define structure-geometry parameters
	set NStories 2;						# number of stories
	set NBays 1;						# number of frame bays (excludes bay for P-delta column)
	set WBay [expr 30.0*12.0];			# bay width in inches
	set HStory1 [expr 15.0*12.0];		# 1st story height in inches
	set HStoryTyp [expr 12.0*12.0];		# story height of other stories in inches
	set HBuilding [expr $HStory1 + ($NStories-1)*$HStoryTyp];	# height of building

# calculate locations of beam/column joints:
	set Pier1  0.0;		# leftmost column line
	set Pier2  [expr $Pier1 + $WBay];
	set Pier3  [expr $Pier2 + $WBay];	# P-delta column line	
	set Floor1 0.0;	# ground floor
	set Floor2 [expr $Floor1 + $HStory1];
	set Floor3 [expr $Floor2 + $HStoryTyp];

# calculate joint offset distance for beam plastic hinges
	set phlat23 [expr 0.0];		# lateral dist from beam-col joint to loc of hinge on Floor 2

# calculate nodal masses -- lump floor masses at frame nodes
	set g 386.2;				# acceleration due to gravity
	set Floor2Weight 535.0;		# weight of Floor 2 in kips
	set Floor3Weight 525.0;		# weight of Floor 3 in kips
	set WBuilding [expr $Floor2Weight + $Floor3Weight];	# total building weight
	set NodalMass2 [expr ($Floor2Weight/$g) / (2.0)];	# mass at each node on Floor 2
	set NodalMass3 [expr ($Floor3Weight/$g) / (2.0)];	# mass at each node on Floor 3
	set Negligible 1e-9;	# a very small number to avoid problems with zero

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
	
# extra nodes for springs in p-delta columns 
	node 326 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	node 327 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	node 336 $Pier3 $Floor3;	# zero-stiffness spring will be used on p-delta column

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
	set Es 29000.0;		# steel Young's modulus

# define column section W24x131 for Story 1 & 2
	set Acol_12 38.5;		# cross-sectional area
	set Icol_12 4020.0;		# moment of inertia
	set Mycol_12 20350.0;	# yield moment

# define beam section W27x102 for Floor 2 & 3
	set Abeam_23 30.0;		# cross-sectional area (full section properties)ce
	set Ibeam_23 3620.0;	# moment of inertia  (full section properties)
	set Mybeam_23 10938.0;	# yield moment at plastic hinge location (i.e., My of RBS section, if used)
	# note: In this example the hinges form right at the beam-column joint, so using an RBS doesn't make sense; 
	#		however, it is done here simply for illustrative purposes.
	
# define plastic hinge lengths
	set Lp_c1 [expr 0.004*$HStory1];		# length of plastic hinge for Story 1 columns
	set Lp_cTyp [expr 0.004*$HStoryTyp];	# length of plastic hinge for typical columns
	set Lp_b23 [expr 0.004*$WBay];			# length of plastic hinge for Floor 2 & 3 beams
		
# determine stiffness modifications so that the strain hardening of the plastic hinge region captures the actual frame member's strain hardening
	# Reference:  Ibarra, L. F., and Krawinkler, H. (2005). "Global collapse of frame structures under seismic excitations," Technical Report 152,
	#             The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
	set n_c1 [expr $HStory1/$Lp_c1];		# rotational stiffness ratio:  (Story 1 column plastic hinge region) / (actual Story 1 column)
	set n_cTyp [expr $HStoryTyp/$Lp_cTyp];	# rotational stiffness ratio:  (Story 2 column plastic hinge region) / (actual Story 2 column)
	set n_b23 [expr $WBay/$Lp_b23];			# rotational stiffness ratio:  (beam plastic hinge region) / (actual beam)
	
# calculate rotational stiffness for plastic hinges
	set Ks_col_1 [expr 6.0*$Es*$Icol_12/$Lp_c1];	# rotational stiffness of Story 1 column hinges
	set Ks_col_2 [expr 6.0*$Es*$Icol_12/$Lp_cTyp];	# rotational stiffness of Story 2 column hinges
	set Ks_beam_23 [expr 6.0*$Es*$Ibeam_23/$Lp_b23];# rotational stiffness of Floor 2 & 3 beam hinges
	
	set Kmem_col_1 [expr 6.0*$Es*$Icol_12/$HStory1];	# rotational stiffness of Story 1 columns
	set Kmem_col_2 [expr 6.0*$Es*$Icol_12/$HStoryTyp];	# rotational stiffness of Story 2 columns
	set Kmem_beam_23 [expr 6.0*$Es*$Ibeam_23/$WBay];	# rotational stiffness of Floor 2 & 3 beams
	
###################################################################################################
#          Define Rotational Springs for Plastic Hinges												  
###################################################################################################
# define rotational spring properties and create spring elements using "rotSect2DModIKModel" procedure
	# rotSect2DModIKModel creates a section with an elastic axial and bilinear flexural response based on Modified Ibarra Krawinkler Deterioration Model
	# references provided in rotSect2DModIKModel.tcl
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
	set a_mem [expr ($Mycol_12*($McMy-1.0)) / ($Kmem_col_1*$th_pP)];	# strain hardening ratio of member
	set b [expr ($a_mem)/(1.0+$n_c1*(1.0-$a_mem))];						# modified strain hardening ratio (Ibarra & Krawinkler 2005, note: Eqn B.5 is incorrect)

# define column plastic hinge sections
	# command: rotSect2DModIKModel	id    ndR  ndC     K   asPos  asNeg  MyPos      MyNeg      LS    LK    LA    LD   cS   cK   cA   cD  th_p+   th_p-   th_pc+   th_pc-  Res+   Res-   th_u+   th_u-    D+     D-
	# Section ID: "10y", where 10 = col section, y = Story #
	# Story 1 columns
	set sec_c1 101;		# section ID for Story 1 columns
	rotSect2DModIKModel $sec_c1 $Es $Acol_12 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
 
	# recompute strain hardening since Story 2 is not the same height as Story 1
	set a_mem [expr ($Mycol_12*($McMy-1.0)) / ($Kmem_col_2*$th_pP)];	# strain hardening ratio of member
	set b [expr ($a_mem)/(1.0+$n_cTyp*(1.0-$a_mem))];       			# modified strain hardening ratio (Ibarra & Krawinkler 2005, note: there is mistake in Eqn B.5)
	# Story 2 columns
	set sec_cTyp 102;	# section ID for Story 2 columns
	rotSect2DModIKModel $sec_cTyp $Es $Acol_12 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	
# define beam plastic hinge sections
	# Setion ID: "20y", where 20 = beam section, y = Floor #
	# redefine the rotations since they are not the same
	set th_pP 0.02;
	set th_pN 0.02;
	set th_pcP 0.16;
	set th_pcN 0.16;
	set a_mem [expr ($Mybeam_23*($McMy-1.0)) / ($Kmem_beam_23*$th_pP)];	# strain hardening ratio of member
	set b [expr ($a_mem)/(1.0+$n_b23*(1.0-$a_mem))];       				# modified strain hardening ratio (Ibarra & Krawinkler 2005, note: there is mistake in Eqn B.5)
	#beam sections
	set sec_b23 202;	# section ID for beams
	rotSect2DModIKModel $sec_b23 $Es $Abeam_23 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

# set up geometric transformations of element
	set PDeltaTransf 1;
	geomTransf PDelta $PDeltaTransf; 	# PDelta transformation

# define column elements using "element" command
	# command: element beamWithHinges $eleID $iNode $jNode $secTagI $Lpi $secTagJ $Lpj $E $A $I $transfID
	# eleID convention:  "1xy" where 1 = col, x = Pier #, y = Story #
	# Columns Story 1
	element beamWithHinges  111  11 12 $sec_c1 $Lp_c1 $sec_c1 $Lp_c1 $Es $Acol_12  $Icol_12 $PDeltaTransf;	# Pier 1
	element beamWithHinges  121  21 22 $sec_c1 $Lp_c1 $sec_c1 $Lp_c1 $Es $Acol_12  $Icol_12 $PDeltaTransf;	# Pier 2
	# Columns Story 2
	element beamWithHinges  112  12 13 $sec_cTyp $Lp_cTyp $sec_cTyp $Lp_cTyp $Es $Acol_12 $Icol_12 $PDeltaTransf;	# Pier 1
	element beamWithHinges  122  22 23 $sec_cTyp $Lp_cTyp $sec_cTyp $Lp_cTyp $Es $Acol_12 $Icol_12 $PDeltaTransf;	# Pier 2
	
# define elastic beam elements
	# eleID convention:  "2xy" where 2 = beam, x = Bay #, y = Floor #
	# Beams Story 1
	element beamWithHinges  212  12 22 $sec_b23 $Lp_b23 $sec_b23 $Lp_b23  $Es $Abeam_23 $Ibeam_23 $PDeltaTransf;
	# Beams Story 2
	element beamWithHinges  222  13 23 $sec_b23 $Lp_b23 $sec_b23 $Lp_b23  $Es $Abeam_23 $Ibeam_23 $PDeltaTransf;
	
# define p-delta columns and rigid links
	set TrussMatID 600;		# define a material ID
	set Arigid 1000.0;		# define area of truss section (make much larger than A of frame elements)
	set Irigid 100000.0; 	# moment of inertia for p-delta columns  (make much larger than I of frame elements)
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
	
# define p-delta column spring: zero-stiffness elastic spring	
	#Spring ID: "5xya" where 5 = leaning column spring, x = Pier #, y = Story #, a = location in story
	# "a" convention: 1 = bottom of story, 2 = top of story
	# rotLeaningCol	ElemID ndR ndC 
	rotLeaningCol 5312 32 326;	# top of Story 1
	rotLeaningCol 5321 32 327;	# bottom of Story 2
	rotLeaningCol 5322 33 336;	# top of Story 2
	
# display the model with the node numbers
	DisplayModel2D NodeNumbers
	
############################################################################
#                       Eigenvalue Analysis                    			   
############################################################################
	set pi [expr 2.0*asin(1.0)];              			# Definition of pi
	set nEigenI 1;										# mode i = 1
	set nEigenJ 2;										# mode j = 2
	set lambdaN [eigen [expr $nEigenJ]];				# eigenvalue analysis for nEigenJ modes
	set lambdaI [lindex $lambdaN [expr 0]]; 			# eigenvalue mode i = 1
	set lambdaJ [lindex $lambdaN [expr $nEigenJ-1]]; 	# eigenvalue mode j = 2
	set w1 [expr pow($lambdaI,0.5)];					# w1 (1st mode circular frequency)
	set w2 [expr pow($lambdaJ,0.5)];					# w2 (2nd mode circular frequency)
	set T1 [expr 2.0*$pi/$w1];      					# 1st mode period of the structure
	set T2 [expr 2.0*$pi/$w2];      					# 2nd mode period of the structure
	puts "T1 = $T1 s"
	puts "T2 = $T2 s"
	
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
	variable constraintsTypeGravity Plain;	# default;
	constraints $constraintsTypeGravity ;   # how it handles boundary conditions
	numberer RCM;							# renumber dof's to minimize band-width (optimization)
	system BandGeneral ;					# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	test NormDispIncr $Tol 6 ; 				# determine if convergence has been achieved at the end of an iteration step
	algorithm Newton;						# use Newton's solution algorithm: updates tangent stiffness at every iteration
	set NstepGravity 10;  					# apply gravity in 10 steps
	set DGravity [expr 1.0/$NstepGravity]; 	# first load increment;
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
	recorder Node -file $dataDir/Vbase.out -node 11 21 31 -dof 1 reaction;	
	
# record story 1 column forces in global coordinates 
	recorder Element -file $dataDir/Fcol111.out -ele 111 force;
	recorder Element -file $dataDir/Fcol121.out -ele 121 force;
	recorder Element -file $dataDir/Fcol731.out -ele 731 force;

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