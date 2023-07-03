# --------------------------------------------------------------------------------------------------
# Example: 1-Story 1-Bay Special Concentric Braced Frame 
# Braces are defined with two force beam-column elements (with 5 IPs) 
# Elements are fully rigid when in thuch with gusset plates
# Rigid liks of the braces are defined with rigid elastic elements. Beam and column reigid elemnts are defined using rigid offsets inside of geometric transformation. 

# Created by:  Vesna Terzic, UC Berkeley, 2011
# Units: kips, inches, seconds

# Element and Node ID conventions:
#	1xy = frame columns 
#	2xy = frame beams 
#	3xy = braces
#	6xy = rigid links

		
###################################################################################################
#          Set Up & Source Definition									  
###################################################################################################
	wipe;							    # clear memory of past model definitions
	model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm = #dimension, ndf = #dofs
	source WSection.tcl;				# procedure for creating standard steel W section
	source HSSsection.tcl;				# procedure for creating standard steel HSS section
	
###################################################################################################
#          Define Building Geometry, Nodes, and Constraints											  
###################################################################################################
# define structure-geometry parameters
	set NStories 1;						# number of stories
	set NBays 1;						# number of frame bays (excludes bay for P-delta column)
	set WBay      360.;		            # bay width in inches
	set HStory    180.;         		# 1st story height in inches
	set HBuilding [expr $NStories*$HStory];	# height of building

# calculate locations frame nodes:
	set Pier1  0.0;						# leftmost column line
	set Pier2  [expr $Pier1 + $WBay];
	set Floor1 0.0;						# ground floor
	set Floor2 [expr $Floor1 + $HStory];


# joint offset distance for beams columns and braces (it is assumed that the the rigid brace lenght is 0.15*Lbr and the angle alpha for a gusset plate is 25 degrees)
	set jOff_col_1 27.18;    	# joint offset for columns at the Storey 1
	set jOff_beam_2 29.25;    	# joint offset for beams at Floor 2
	set jOff_braceX  27.0;    	# joint offset for all braces in X direction
	set jOff_braceY  27.0;    	# joint offset for all braces in Y direction

# calculate nodal masses -- lump floor masses at frame nodes
	set g 386.4;			# acceleration due to gravity
	# horizontal mass (half of the floor mass devided on two columns)
	set NodalMass2H 1.26;	# mass at each column node on Floor 2
	# veritcal mass (tributary mass around the column: http://onlinelibrary.wiley.com/doi/10.1002/tal.251/pdf)
	set NodalMass2V 0.105;  # mass at each column node on Floor 2

# define nodes and assign masses to beam-column intersections of frame
	# command:  node nodeID xcoord ycoord -mass mass_dof1 mass_dof2 mass_dof3
	# nodeID convention:  "xy" where x = Pier # and y = Floor # 
	node 11 $Pier1 $Floor1;
	node 21 $Pier2 $Floor1;
	node 12 $Pier1 $Floor2 -mass $NodalMass2H $NodalMass2V 0.0;
	node 22 $Pier2 $Floor2 -mass $NodalMass2H $NodalMass2V 0.0;	 

# define extra nodes in the braces
	# nodeID convention:  "1xa" where x = storey # and a = location of the node
	# "a" convention: 1 = left; 2 = right;
	# nodes at Storey 1 
	node 111 [expr $Pier1+$WBay/4.0] [expr $Floor1+$HStory/2.0]
	node 112 [expr $Pier2-$WBay/4.0] [expr $Floor1+$HStory/2.0]
	
# define extra nodes in the beams
	# nodeID convention:  "2xa" where x = floor # and a = location of the node
	# "a" convention: 1 = left; 2 = middle; 3 = right;
	# nodes at Floor 2 
	node 221 $Pier1 $Floor2
	node 222 [expr $WBay/2.0] $Floor2
	node 223 $Pier2 $Floor2
	
# define extra nodes for rigid links in the braces:
	# nodeID convention:  "4xa" where x = storey # and a = location of the node
	# "a" convention: 1 = left-botom; 2 = left-top; 3 = right-bottom; 4 = right-top;
	# nodes at storey 1
	node 411 [expr $Pier1+$jOff_braceX] 			[expr $Floor1+$jOff_braceY]
	node 412 [expr $Pier1+$WBay/2.0-$jOff_braceX] 	[expr $Floor2-$jOff_braceY]
	node 413 [expr $Pier2-$jOff_braceX] 			[expr $Floor1+$jOff_braceY]
	node 414 [expr $Pier2-$WBay/2.0+$jOff_braceX] 	[expr $Floor2-$jOff_braceY]  

# assign boundary condidtions 
	# command:  fix nodeID dxFixity dyFixity rzFixity
	# fixity values: 1 = constrained; 0 = unconstrained
	# fix the base of the building;
	fix 11 1 1 1;
	fix 21 1 1 1;
	
# define constraints for pined beam-to-column connection
	# beams of the floor 2
	equalDOF 12 221 1 2 
	equalDOF 22 223 1 2 

###################################################################################################
#          Define Materials and Sections													  
###################################################################################################

# define material for nonlinear beams and columns
	set matID_BC 1
	set matID_fatBC 2
	set Es 29000.0;  # modulus of elasticity for steel
	set Fy 50.0; 	 # yield stress of steel
	set b 0.003;	 # strain hardening ratio
	uniaxialMaterial Steel02 $matID_BC $Fy $Es $b 20 0.925 0.15
	uniaxialMaterial Fatigue $matID_fatBC $matID_BC 


# define material for braces 
	set matID_Brace 3
	set matID_fatBrace 4
	set E0 0.095
	set m -0.5
	uniaxialMaterial Steel02 $matID_Brace $Fy $Es $b 20 0.925 0.15 0.0005 0.01 0.0005 0.01
	uniaxialMaterial Fatigue $matID_fatBrace $matID_Brace -E0 $E0 -m $m

# define sections for columns and beams
	# secID: "ax", x=storey # or floor #, a = type of element
	# "a" convention: 1 = column; 2 = beam;
	# define column sections
	# command: WSection  secID matID d bf tf tw nfdw nftw nfbf nftf
	# W14x176 
	WSection 11 $matID_fatBC 15.22 15.65  1.31 0.83  14 2 14 2 
	# define beam sections
 	# W27x84
	WSection 22 $matID_fatBC 26.71 9.96   0.64 0.46  20 2 10 2 
	
# define sections for braces
	# secID convention: "3x", x=storey #
	# command: HSSsection secID matID d t nfdy nfty nfdz nftz
	# HSS12x12x0.625
	HSSsection 31 $matID_fatBrace 12. 0.625 12 2 12 2  
		
###################################################################################################
#          Define Geometric Transformation													  
###################################################################################################

# columns
	# columns of a braced frame
	# transfTag convention: "1x", x=storey #
	geomTransf PDelta 11  -jntOffset 0.0 $jOff_col_1 0.0 0.0

# beams
	# transfTag convention: "2xa", x=Floor #, a = location of the beam
	# "a" convention: 1 = left; 2 = right;
	geomTransf PDelta 221  -jntOffset 0.0 0.0 -$jOff_beam_2 0.0
	geomTransf PDelta 222  -jntOffset $jOff_beam_2 0.0 0.0 0.0
	
# braces
	geomTransf Corotational 31

# rigid links
	geomTransf Linear 41

	puts "A"
###################################################################################################
#          Define Elements													  
###################################################################################################

# define columns of a braced frame:
	# eleID convention: "1xy", 1 = column, x=Pier #, y= storey #
	# command arguments:     $eleID $iNode $jNode $numIntgrPts $secTag $transfTag
	# storey 1	
	element forceBeamColumn   111      11     12        5         11       11
	element forceBeamColumn   121      21     22        5         11       11	
puts "B"	
# define beams of a braced frame:
	# eleID convention: "2xa", 2 = beam, x=floor #, a = location of the beam
	# "a" convention: 1 = left; 2 = right; 0 = the whole beam
	# command arguments:     $eleID $iNode $jNode $numIntgrPts $secTag $transfTag
	# floor 2
	element forceBeamColumn   221     221    222       5         22       221
	element forceBeamColumn   222     222    223       5         22       222
puts "C"
# define braces:
	# eleID convention: "3xa", 3 = brace, x=storey #, a = location of the beam
	# "a" convention: 1 = the 1st from left; 2 = the 2nd from left; 3 = the 1st from right; 4 = the 1st from right;
	# command arguments:     $eleID $iNode $jNode $numIntgrPts $secTag $transfTag
	# storey 1	
	element forceBeamColumn   311     411    111        5         31       31
	element forceBeamColumn   312     111    412        5         31       31
	element forceBeamColumn   313     413    112        5         31       31
	element forceBeamColumn   314     112    414        5         31       31
puts "D"
	
# define rigid links:
# 	set Arigid1 25.7
# 	set Arigid2 21.
# 	set Arigid3 13.5
	set Arigid1 [expr 25.7*10.]
	set Arigid2 [expr 21.*10.]
	set Arigid3 [expr 13.5*10.]
	set Irigid 30000.
	# eleID convention: "6xa", 6 = rigid link, x = storey #, a = location of the rigid link
	# "a" convention: 1 = left-botom; 2 = left-top; 3 = right-bottom; 4 = right-top;
	# comand arguemnts:       $eleTag $iNode $jNode   $A     $E    $Iz   $transfTag
	# storey 1
	element elasticBeamColumn  	611     11     411 	$Arigid1  $Es  $Irigid     41
	element elasticBeamColumn  	612    222     412 	$Arigid1  $Es  $Irigid     41
	element elasticBeamColumn  	613     21     413 	$Arigid1  $Es  $Irigid     41
	element elasticBeamColumn  	614    222     414 	$Arigid1  $Es  $Irigid     41
	
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

	recorder display "Mode 1" 10 10 500 500 -wipe 
	prp 360. 270. 1.;
	vup  0  1 0;
	vpn  0  0 1;
	display -1 2 50
 		
############################################################################
#              Gravity Loads & Gravity Analysis
############################################################################
# apply gravity loads
	#command: pattern PatternType $PatternID TimeSeriesType

	pattern Plain 101 Linear {
		
		# point loads on frame column nodes
		set P_F2    -40.46;	# load on each frame node in Floor 2 

		# Floor 2 loads
		load 12 0.0 $P_F2 0.0;
		load 22 0.0 $P_F2 0.0;
	}

# Gravity-analysis: load-controlled static analysis
	set Tol 1.0e-6;							# convergence tolerance for test
	constraints Transformation;				# how it handles boundary conditions
	#constraints Penalty 1.0e15 1.0e15;
	numberer RCM;							# renumber dof's to minimize band-width (optimization)
	system BandGeneral;						# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	test NormDispIncr $Tol 10;				# determine if convergence has been achieved at the end of an iteration step
	algorithm Newton ;						# use Newton's solution algorithm: updates tangent stiffness at every iteration
	set NstepGravity 10;					# apply gravity in 10 steps
	set DGravity [expr 1.0/$NstepGravity];	# load increment
	integrator LoadControl $DGravity;		# determine the next time step for an analysis
	analysis Static;						# define type of analysis static or transient
	analyze $NstepGravity;					# apply gravity

	# maintain constant gravity loads and reset time to zero
	loadConst -time 0.0
	puts "Model Built"
	
