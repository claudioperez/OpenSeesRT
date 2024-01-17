# --------------------------------------------------------------------------------------------------
# Example: 2-Story Steel Moment Frame with Concentrated Plasticity
# Panel Zone Model with Concentrated Plastic Hinges at RBS Locations
# Created by:  Laura Eads, Stanford University, 2010
# Units: kips, inches, seconds

# Updated 9 May 2012:  fixed errors defining rayleigh damping (see line 556)

# Element ID conventions:
#	1xy    = frame columns with RBS springs at both ends
#	2xy    = frame beams with RBS springs at both ends
#	6xy    = trusses linking frame and P-delta column
#	7xy    = P-delta columns
#	2xya   = frame beams between panel zone and RBS spring
#	3xya   = frame column rotational springs
#	4xya   = frame beam rotational springs
#	5xya   = P-delta column rotational springs
#	4xy00  = panel zone rotational springs
#	500xya = panel zone elements
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
	source RotSpring2DModIKModel.tcl;	# procedure for defining a rotational spring (zero-length element) for plastic hinges
	source RotLeaningCol.tcl;			# procedure for defining a rotational spring (zero-length element) with very small stiffness
	source RotPanelZone2D.tcl;			# procedure for defining a rotational spring (zero-length element) to capture panel zone shear distortions
	source ElemPanelZone2D.tcl;			# procedure for defining 8 elements to create a rectangular panel zone

###################################################################################################
#          Define Analysis Type
###################################################################################################
# Define type of analysis:  "pushover" = pushover;  "dynamic" = dynamic
	set analysisType "dynamic";

	if {$analysisType == "pushover"} {
		set dataDir Concentrated-PanelZone-Pushover-Output;	# name of output folder
	}
	if {$analysisType == "dynamic"} {
		set dataDir Concentrated-PanelZone-Dynamic-Output;	# name of output folder
	}
	file mkdir $dataDir; 					        # create output folder

###################################################################################################
#          Define Building Geometry, Nodes, Masses, and Constraints
###################################################################################################
# define structure-geometry parameters
	set NStories 2;						# number of stories
	set NBays 1;						# number of frame bays (excludes bay for P-delta column)
	set WBay      [expr 30.0*12.0];		# bay width in inches
	set HStory1   [expr 15.0*12.0];		# 1st story height in inches
	set HStoryTyp [expr 12.0*12.0];		# story height of other stories in inches
	set HBuilding [expr $HStory1 + ($NStories-1)*$HStoryTyp];	# height of building

# calculate locations of beam-column joint centerlines:
	set Pier1  0.0;		# leftmost column line
	set Pier2  [expr $Pier1 + $WBay];
	set Pier3  [expr $Pier2 + $WBay];	# P-delta column line
	set Floor1 0.0;		# ground floor
	set Floor2 [expr $Floor1 + $HStory1];
	set Floor3 [expr $Floor2 + $HStoryTyp];

# calculate panel zone dimensions
	set pzlat23  [expr 24.5/2.0];	# lateral dist from CL of beam-col joint to edge of panel zone (= half the column depth)
	set pzvert23 [expr 27.1/2.0];	# vert dist from CL of beam-col joint to edge of panel zone (= half the beam depth)

# calculate plastic hinge offsets from beam-column centerlines:
	set phlat23 [expr $pzlat23 + 7.5 + 22.5/2.0];	# lateral dist from CL of beam-col joint to beam hinge
	set phvert23 [expr $pzvert23 + 0.0];			# vert dist from CL of beam-col joint to column hinge (forms at edge of panel zone)

# calculate nodal masses -- lump floor masses at frame nodes
	set g 386.2;				# acceleration due to gravity
	set Floor2Weight 535.0;		# weight of Floor 2 in kips
	set Floor3Weight 525.0;		# weight of Floor 3 in kips
	set WBuilding  [expr $Floor2Weight + $Floor3Weight];# total building weight
	set NodalMass2 [expr ($Floor2Weight/$g) / (2.0)];	# mass at each node on Floor 2
	set NodalMass3 [expr ($Floor3Weight/$g) / (2.0)];	# mass at each node on Floor 3
	set Negligible 1e-9;	# a very small number to avoid problems with zero

# define nodes and assign masses to beam-column intersections of frame
	# command:  node nodeID xcoord ycoord -mass mass_dof1 mass_dof2 mass_dof3
	# nodeID convention:  "xy" where x = Pier # and y = Floor #
	node 11 $Pier1 $Floor1;
	node 21 $Pier2 $Floor1;
	node 31 $Pier3 $Floor1;
	node 32 $Pier3 $Floor2;
	node 33 $Pier3 $Floor3;

# define extra nodes for plastic hinge rotational springs
	# nodeID convention:  "xya" where x = Pier #, y = Floor #, a = location relative to beam-column joint
	# "a" convention: 1,2 = left; 3,4 = right; (used for beams)
	# "a" convention: 5,6 = below; 7,8 = above; (used for columns)
	# column hinges at bottom of Story 1 (base)
	node 117 $Pier1 $Floor1;
	node 217 $Pier2 $Floor1;
	# column hinges at top of Story 1
	node 125 $Pier1 [expr $Floor2 - $phvert23];
	node 126 $Pier1 [expr $Floor2 - $phvert23];
	node 225 $Pier2 [expr $Floor2 - $phvert23];
	node 226 $Pier2 [expr $Floor2 - $phvert23];
	node 326 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	# column hinges at bottom of Story 2
	node 127 $Pier1 [expr $Floor2 + $phvert23];
	node 128 $Pier1 [expr $Floor2 + $phvert23];
	node 227 $Pier2 [expr $Floor2 + $phvert23];
	node 228 $Pier2 [expr $Floor2 + $phvert23];
	node 327 $Pier3 $Floor2;	# zero-stiffness spring will be used on p-delta column
	# column hinges at top of Story 2
	node 135 $Pier1 [expr $Floor3 - $phvert23];
	node 136 $Pier1 [expr $Floor3 - $phvert23];
	node 235 $Pier2 [expr $Floor3 - $phvert23];
	node 236 $Pier2 [expr $Floor3 - $phvert23];
	node 336 $Pier3 $Floor3;	# zero-stiffness spring will be used on p-delta column

	# beam hinges at Floor 2
	node 121 [expr $Pier1 + $phlat23] $Floor2;
	node 122 [expr $Pier1 + $phlat23] $Floor2;
	node 223 [expr $Pier2 - $phlat23] $Floor2;
	node 224 [expr $Pier2 - $phlat23] $Floor2;
	# beam hinges at Floor 3
	node 131 [expr $Pier1 + $phlat23] $Floor3;
	node 132 [expr $Pier1 + $phlat23] $Floor3;
	node 233 [expr $Pier2 - $phlat23] $Floor3;
	node 234 [expr $Pier2 - $phlat23] $Floor3;

# define extra nodes for panel zones
	# nodeID convention:  "xybc" where x = Pier #, y = Floor #, bc = location relative to beam-column joint
	# "bc" conventions: 01,02 = top left of joint;
	# 					03,04 = top right of joint;
	# 					05    = middle right of joint; (vertical middle, horizontal right)
	# 					06,07 = btm right of joint;
	# 					08,09 = btm left of joint;
	# 					10    = middle left of joint; (vertical middle, horizontal left)
	# note: top center and btm center nodes were previously defined as xy7 and xy6, respectively, at Floor 2(center = horizonal center)

	# panel zone at Pier 1, Floor 2
	node 1201 [expr $Pier1 - $pzlat23] [expr $Floor2 + $phvert23];
	node 1202 [expr $Pier1 - $pzlat23] [expr $Floor2 + $phvert23];
	node 1203 [expr $Pier1 + $pzlat23] [expr $Floor2 + $phvert23];
	node 1204 [expr $Pier1 + $pzlat23] [expr $Floor2 + $phvert23];
	node 1205 [expr $Pier1 + $pzlat23] [expr $Floor2];
	node 1206 [expr $Pier1 + $pzlat23] [expr $Floor2 - $phvert23];
	node 1207 [expr $Pier1 + $pzlat23] [expr $Floor2 - $phvert23];
	node 1208 [expr $Pier1 - $pzlat23] [expr $Floor2 - $phvert23];
	node 1209 [expr $Pier1 - $pzlat23] [expr $Floor2 - $phvert23];
	node 1210 [expr $Pier1 - $pzlat23] [expr $Floor2];

	# panel zone at Pier 2, Floor 2
	node 2201 [expr $Pier2 - $pzlat23] [expr $Floor2 + $phvert23];
	node 2202 [expr $Pier2 - $pzlat23] [expr $Floor2 + $phvert23];
	node 2203 [expr $Pier2 + $pzlat23] [expr $Floor2 + $phvert23];
	node 2204 [expr $Pier2 + $pzlat23] [expr $Floor2 + $phvert23];
	node 2205 [expr $Pier2 + $pzlat23] [expr $Floor2];
	node 2206 [expr $Pier2 + $pzlat23] [expr $Floor2 - $phvert23];
	node 2207 [expr $Pier2 + $pzlat23] [expr $Floor2 - $phvert23];
	node 2208 [expr $Pier2 - $pzlat23] [expr $Floor2 - $phvert23];
	node 2209 [expr $Pier2 - $pzlat23] [expr $Floor2 - $phvert23];
	node 2210 [expr $Pier2 - $pzlat23] [expr $Floor2];

	# panel zone at Pier 1, Floor 3
	node 1301 [expr $Pier1 - $pzlat23] [expr $Floor3 + $phvert23];
	node 1302 [expr $Pier1 - $pzlat23] [expr $Floor3 + $phvert23];
	node 1303 [expr $Pier1 + $pzlat23] [expr $Floor3 + $phvert23];
	node 1304 [expr $Pier1 + $pzlat23] [expr $Floor3 + $phvert23];
	node 1305 [expr $Pier1 + $pzlat23] [expr $Floor3];
	node 1306 [expr $Pier1 + $pzlat23] [expr $Floor3 - $phvert23];
	node 1307 [expr $Pier1 + $pzlat23] [expr $Floor3 - $phvert23];
	node 1308 [expr $Pier1 - $pzlat23] [expr $Floor3 - $phvert23];
	node 1309 [expr $Pier1 - $pzlat23] [expr $Floor3 - $phvert23];
	node 1310 [expr $Pier1 - $pzlat23] [expr $Floor3];
	node 137  [expr $Pier1]  [expr $Floor3 + $phvert23]; # not previously defined since no column above

	# panel zone at Pier 2, Floor 3
	node 2301 [expr $Pier2 - $pzlat23] [expr $Floor3 + $phvert23];
	node 2302 [expr $Pier2 - $pzlat23] [expr $Floor3 + $phvert23];
	node 2303 [expr $Pier2 + $pzlat23] [expr $Floor3 + $phvert23];
	node 2304 [expr $Pier2 + $pzlat23] [expr $Floor3 + $phvert23];
	node 2305 [expr $Pier2 + $pzlat23] [expr $Floor3];
	node 2306 [expr $Pier2 + $pzlat23] [expr $Floor3 - $phvert23];
	node 2307 [expr $Pier2 + $pzlat23] [expr $Floor3 - $phvert23];
	node 2308 [expr $Pier2 - $pzlat23] [expr $Floor3 - $phvert23];
	node 2309 [expr $Pier2 - $pzlat23] [expr $Floor3 - $phvert23];
	node 2310 [expr $Pier2 - $pzlat23] [expr $Floor3];
	node 237  [expr $Pier2]  [expr $Floor3 + $phvert23]; # not previously defined since no column above

# define nodal masses:  lump at beam-column joints in frame
	# command: mass $nodeID $dof1mass $dof2mass $dof3mass
	mass 1205 $NodalMass2 $Negligible $Negligible;	# Pier 1, Floor 2
	mass 2205 $NodalMass2 $Negligible $Negligible;	# Pier 2, Floor 2
	mass 1305 $NodalMass3 $Negligible $Negligible;	# Pier 1, Floor 3
	mass 2305 $NodalMass3 $Negligible $Negligible;	# Pier 2, Floor 3

# constrain beam-column joints in a floor to have the same lateral displacement using the "equalDOF" command
	# command: equalDOF $MasterNodeID $SlaveNodeID $dof1 $dof2...
	set dof1 1;	# constrain movement in dof 1 (x-direction)
	equalDOF 1205 2205 $dof1;	# Floor 2:  Pier 1 to Pier 2
	equalDOF 1205 32 $dof1;		# Floor 2:  Pier 1 to Pier 3
	equalDOF 1305 2305 $dof1;	# Floor 3:  Pier 1 to Pier 2
	equalDOF 1305 33 $dof1;		# Floor 3:  Pier 1 to Pier 3

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
	set Fy 50.0;			# steel yield strength

# define column section W24x131 for Story 1 & 2
	set Acol_12  38.5;		# cross-sectional area
	set Icol_12  4020.0;	# moment of inertia
	set Mycol_12 20350.0;	# yield moment at plastic hinge location (i.e., My of RBS section)
	set dcol_12 24.5;		# depth
	set bfcol_12 12.9;		# flange width
	set tfcol_12 0.96;		# flange thickness
	set twcol_12 0.605;		# web thickness

# define beam section W27x102 for Floor 2 & 3
	set Abeam_23  30.0;		# cross-sectional area (full section properties)
	set Ibeam_23  3620.0;	# moment of inertia  (full section properties)
	set Mybeam_23 10938.0;	# yield moment at plastic hinge location (i.e., My of RBS section)
	set dbeam_23 27.1;		# depth

# determine stiffness modifications to equate the stiffness of the spring-elastic element-spring subassembly to the stiffness of the actual frame member
	# References: (1) Ibarra, L. F., and Krawinkler, H. (2005). "Global collapse of frame structures under seismic excitations," Technical Report 152,
	#             		The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
	#			  (2) Zareian, F. and Medina, R. A. (2010). “A practical method for proper modeling of structural damping in inelastic plane
	#					structural systems,” Computers & Structures, Vol. 88, 1-2, pp. 45-53.
	# calculate modified section properties to account for spring stiffness being in series with the elastic element stiffness
	set n 10.0;		# stiffness multiplier for rotational spring

	# calculate modified moment of inertia for elastic elements between plastic hinge springs
	set Icol_12mod  [expr $Icol_12*($n+1.0)/$n];	# modified moment of inertia for columns in Story 1 & 2
	set Ibeam_23mod [expr $Ibeam_23*($n+1.0)/$n];	# modified moment of inertia for beams in Floor 2 & 3

	# calculate modified rotational stiffness for plastic hinge springs: use length between springs
	set Ks_col_1   [expr $n*6.0*$Es*$Icol_12mod/($HStory1-$phvert23)];		# rotational stiffness of Story 1 column springs
	set Ks_col_2   [expr $n*6.0*$Es*$Icol_12mod/($HStoryTyp-2*$phvert23)];	# rotational stiffness of Story 2 column springs
	set Ks_beam_23 [expr $n*6.0*$Es*$Ibeam_23mod/($WBay-2*$phlat23)];		# rotational stiffness of Floor 2 & 3 beam springs

# set up geometric transformation of elements
	set PDeltaTransf 1;
	geomTransf PDelta $PDeltaTransf; 	# PDelta transformation

# define elastic column elements using "element" command
	# command: element elasticBeamColumn $eleID $iNode $jNode $A $E $I $transfID
	# eleID convention:  "1xy" where 1 = col, x = Pier #, y = Story #
	# Columns Story 1
	element elasticBeamColumn  111  117 125 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 1
	element elasticBeamColumn  121  217 225 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 2
	# Columns Story 2
	element elasticBeamColumn  112  128 135 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 1
	element elasticBeamColumn  122  228 235 $Acol_12 $Es $Icol_12mod $PDeltaTransf;	# Pier 2

# define elastic beam elements
	# element between plastic hinges: eleID convention = "2xy" where 2 = beam, x = Bay #, y = Floor #
	# element between plastic hinge and panel zone: eleID convention = "2xya" where 2 = beam, x = Bay #, y = Floor #, a = loc in bay
	#	"a" convention: 1 = left end of bay; 2 = right end of bay
	# Beams Story 1
	element elasticBeamColumn  2121 1205 121  $Abeam_23 $Es $Ibeam_23    $PDeltaTransf;
	element elasticBeamColumn  212  122  223  $Abeam_23 $Es $Ibeam_23mod $PDeltaTransf;
	element elasticBeamColumn  2122 224  2210 $Abeam_23 $Es $Ibeam_23    $PDeltaTransf;
	# Beams Story 2
	element elasticBeamColumn  2131 1305 131  $Abeam_23 $Es $Ibeam_23    $PDeltaTransf;
	element elasticBeamColumn  213  132  233  $Abeam_23 $Es $Ibeam_23mod $PDeltaTransf;
	element elasticBeamColumn  2132 234  2310 $Abeam_23 $Es $Ibeam_23    $PDeltaTransf;

# define p-delta columns and rigid links
	set TrussMatID 600;		# define a material ID
	set Arigid 1000.0;		# define area of truss section (make much larger than A of frame elements)
	set Irigid 100000.0;	# moment of inertia for p-delta columns  (make much larger than I of frame elements)
	uniaxialMaterial Elastic $TrussMatID $Es;		# define truss material
	# rigid links
	# command: element truss $eleID $iNode $jNode $A $materialID
	# eleID convention:  6xy, 6 = truss link, x = Bay #, y = Floor #
	element truss  622 2205 32 $Arigid $TrussMatID;	# Floor 2
	element truss  623 2305 33 $Arigid $TrussMatID;	# Floor 3

	# p-delta columns
	# eleID convention:  7xy, 7 = p-delta columns, x = Pier #, y = Story #
	element elasticBeamColumn  731  31  326 $Arigid $Es $Irigid $PDeltaTransf;	# Story 1
	element elasticBeamColumn  732  327 336 $Arigid $Es $Irigid $PDeltaTransf;	# Story 2

# define elastic panel zone elements (assume rigid)
	# elemPanelZone2D creates 8 elastic elements that form a rectangular panel zone
	# references provided in elemPanelZone2D.tcl
	# note: the nodeID and eleID of the upper left corner of the PZ must be imported
	# eleID convention:  500xya, 500 = panel zone element, x = Pier #, y = Floor #
	# "a" convention: defined in elemPanelZone2D.tcl, but 1 = top left element
	set Apz 1000.0;	# area of panel zone element (make much larger than A of frame elements)
	set Ipz 1.0e5;  # moment of intertia of panel zone element (make much larger than I of frame elements)
	# elemPanelZone2D eleID  nodeR E  A_PZ I_PZ transfTag
	elemPanelZone2D   500121 1201 $Es $Apz $Ipz $PDeltaTransf;	# Pier 1, Floor 2
	elemPanelZone2D   500221 2201 $Es $Apz $Ipz $PDeltaTransf;	# Pier 2, Floor 2
	elemPanelZone2D   500131 1301 $Es $Apz $Ipz $PDeltaTransf;	# Pier 1, Floor 3
	elemPanelZone2D   500231 2301 $Es $Apz $Ipz $PDeltaTransf;	# Pier 2, Floor 3

# display the model with the node numbers
	DisplayModel2D NodeNumbers;

###################################################################################################
#          Define Rotational Springs for Plastic Hinges, Panel Zones, and Leaning Columns
###################################################################################################
# define rotational spring properties and create spring elements using "RotSpring2DModIKModel" procedure
	# RotSpring2DModIKModel creates a uniaxial material spring with a bilinear response based on Modified Ibarra Krawinkler Deterioration Model
	# references provided in RotSpring2DModIKModel.tcl
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
	# command: RotSpring2DModIKModel	id    ndR  ndC     K   asPos  asNeg  MyPos      MyNeg      LS    LK    LA    LD   cS   cK   cA   cD  th_p+   th_p-   th_pc+   th_pc-  Res+   Res-   th_u+   th_u-    D+     D-
	# col springs @ bottom of Story 1 (at base)
	RotSpring2DModIKModel 3111 11 117 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 3211 21 217 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#col springs @ top of Story 1 (below Floor 2)
	RotSpring2DModIKModel 3112 126 125 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 3212 226 225 $Ks_col_1 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

	# recompute strain hardening since Story 2 is not the same height as Story 1
	set a_mem [expr ($n+1.0)*($Mycol_12*($McMy-1.0)) / ($Ks_col_2*$th_pP)];	# strain hardening ratio of spring
	set b [expr ($a_mem)/(1.0+$n*(1.0-$a_mem))];							# modified strain hardening ratio of spring (Ibarra & Krawinkler 2005, note: there is mistake in Eqn B.5)
	# col springs @ bottom of Story 2 (above Floor 2)
	RotSpring2DModIKModel 3121 127 128 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 3221 227 228 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#col springs @ top of Story 2 (below Floor 3)
	RotSpring2DModIKModel 3122 136 135 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 3222 236 235 $Ks_col_2 $b $b $Mycol_12 [expr -$Mycol_12] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

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
	RotSpring2DModIKModel 4121 121 122 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 4122 223 224 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	#beam springs at Floor 3
	RotSpring2DModIKModel 4131 131 132 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;
	RotSpring2DModIKModel 4132 233 234 $Ks_beam_23 $b $b $Mybeam_23 [expr -$Mybeam_23] $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

	# create region for beam springs
	region 2 -ele 4121 4122 4131 4132;

#define panel zone springs
	# rotPanelZone2D creates a uniaxial material spring with a trilinear response based on the Krawinkler Model
	#				It also constrains the nodes in the corners of the panel zone.
	# references provided in rotPanelZone2D.tcl
	# note: the upper right corner nodes of the PZ must be imported
	source RotPanelZone2D.tcl
	set Ry 1.2; 	# expected yield strength multiplier
	set as_PZ 0.03; # strain hardening of panel zones
	# Spring ID: "4xy00" where 4 = panel zone spring, x = Pier #, y = Floor #
	#2nd Floor PZ springs
	#             ElemID  ndR  ndC  E   Fy   dc       bf_c        tf_c       tp        db       Ry   as
	rotPanelZone2D 41200 1203 1204 $Es $Fy $dcol_12 $bfcol_12 $tfcol_12 $twcol_12 $dbeam_23 $Ry $as_PZ;
	rotPanelZone2D 42200 2203 2204 $Es $Fy $dcol_12 $bfcol_12 $tfcol_12 $twcol_12 $dbeam_23 $Ry $as_PZ;
	#3rd Floor PZ springs
	rotPanelZone2D 41300 1303 1304 $Es $Fy $dcol_12 $bfcol_12 $tfcol_12 $twcol_12 $dbeam_23 $Ry $as_PZ;
	rotPanelZone2D 42300 2303 2304 $Es $Fy $dcol_12 $bfcol_12 $tfcol_12 $twcol_12 $dbeam_23 $Ry $as_PZ;

# define p-delta column spring: zero-stiffness elastic spring
	#Spring ID: "5xya" where 5 = leaning column spring, x = Pier #, y = Story #, a = location in story
	# "a" convention: 1 = bottom of story, 2 = top of story
	# rotLeaningCol ElemID ndR ndC
	rotLeaningCol 5312 32 326;	# top of Story 1
	rotLeaningCol 5321 32 327;	# bottom of Story 2
	rotLeaningCol 5322 33 336;	# top of Story 2

	# create region for P-Delta column springs
	region 3 -ele 5312 5321 5322;

        constraints Transformation
############################################################################
#                       Eigenvalue Analysis                    
############################################################################
	set pi [expr 2.0*asin(1.0)];						# Definition of pi
	set nEigenI 1;										# mode i = 1
	set nEigenJ 2;										# mode j = 2
	set lambdaN [eigen [expr $nEigenJ]];				# eigenvalue analysis for nEigenJ modes
	set lambdaI [lindex $lambdaN [expr $nEigenI-1]];	# eigenvalue mode i = 1
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
		load 127 0.0 $P_F2 0.0;
		load 227 0.0 $P_F2 0.0;
		# Floor 3 loads
		load 137 0.0 $P_F3 0.0;
		load 237 0.0 $P_F3 0.0;
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
	analysis Static;						# define type of analysis: static or transient
	analyze $NstepGravity;					# apply gravity

	# maintain constant gravity loads and reset time to zero
	loadConst -time 0.0
	puts "Model Built"

############################################################################
#              Recorders					                
############################################################################
# record drift histories
	# record drifts
	recorder Drift -file $dataDir/Drift-Story1.out -time -iNode 11 -jNode 1205 -dof 1 -perpDirn 2;
	recorder Drift -file $dataDir/Drift-Story2.out -time -iNode 1205 -jNode 1305 -dof 1 -perpDirn 2;
	recorder Drift -file $dataDir/Drift-Roof.out -time -iNode 11 -jNode 1305 -dof 1 -perpDirn 2;

# record floor displacements
	recorder Node -file $dataDir/Disp.out -time -node 1205 1305 -dof 1 disp;

# record base shear reactions
	recorder Node -file $dataDir/Vbase.out -time -node 117 217 31 -dof 1 reaction;

# record story 1 column forces in global coordinates
	recorder Element -file $dataDir/Fcol111.out -time -ele 111 force;
	recorder Element -file $dataDir/Fcol121.out -time -ele 121 force;
	recorder Element -file $dataDir/Fcol731.out -time -ele 731 force;

# record response history of all frame column springs (one file for moment, one for rotation)
	recorder Element -file $dataDir/MRFcol-Mom-Hist.out -time -region 1 force;
	recorder Element -file $dataDir/MRFcol-Rot-Hist.out -time -region 1 deformation;

# record response history of all frame beam springs (one file for moment, one for rotation)
	recorder Element -file $dataDir/MRFbeam-Mom-Hist.out -time -region 2 force;
	recorder Element -file $dataDir/MRFbeam-Rot-Hist.out -time -region 2 deformation;


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
	set lat2 16.255;	# force on each beam-column joint in Floor 2
	set lat3 31.636;	# force on each beam-column joint in Floor 3
	pattern Plain 200 Linear {
					load 1205 $lat2 0.0 0.0;
					load 2205 $lat2 0.0 0.0;
					load 1305 $lat3 0.0 0.0;
					load 2305 $lat3 0.0 0.0;
	}

# display deformed shape:
	set ViewScale 5;
	DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model

# displacement parameters
	set IDctrlNode 1305;				# node where disp is read for disp control
	set IDctrlDOF 1;					# degree of freedom read for disp control (1 = x displacement)
	set Dmax [expr 0.1*$HBuilding];		# maximum displacement of pushover: 10% roof drift
	set Dincr [expr 0.01];				# displacement increment

# analysis commands
	constraints Plain;					# how it handles boundary conditions
	numberer RCM;						# renumber dof's to minimize band-width (optimization)
	system BandGeneral;					# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	test NormUnbalance 1.0e-5 400;		# type of convergence criteria with tolerance, max iterations
	algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
	integrator DisplacementControl  $IDctrlNode   $IDctrlDOF $Dincr;	# use displacement-controlled analysis
	analysis Static;					# define type of analysis: static for pushover
	set Nsteps [expr int($Dmax/$Dincr)];# number of pushover analysis steps
	set ok [analyze $Nsteps];			# this will return zero if no convergence problems were encountered
	puts "Pushover complete";			# display this message in the command window
}

############################################################################
#   Time History/Dynamic Analysis               			   			   #
############################################################################
if {$analysisType == "dynamic"} {
	puts "Running dynamic analysis..."
		# display deformed shape:
		set ViewScale 5;	# amplify display of deformed shape
		DisplayModel2D DeformedShape $ViewScale;	# display deformed shape, the scaling factor needs to be adjusted for each model

	# Rayleigh Damping
		# calculate damping parameters
		set zeta 0.02;		# percentage of critical damping
		set a0 [expr $zeta*2.0*$w1*$w2/($w1 + $w2)];	# mass damping coefficient based on first and second modes
		set a1 [expr $zeta*2.0/($w1 + $w2)];			# stiffness damping coefficient based on first and second modes
		set a1_mod [expr $a1*(1.0+$n)/$n];				# modified stiffness damping coefficient used for n modified elements. See Zareian & Medina 2010.

		# assign damping to frame beams and columns
		# command: region $regionID -eleRange $elementIDfirst $elementIDlast -rayleigh $alpha_mass $alpha_currentStiff $alpha_initialStiff $alpha_committedStiff
		## old commands:	region 4 -eleRange 111 213 rayleigh 0.0 0.0 $a1_mod 0.0;	# assign stiffness proportional damping to frame beams & columns w/ n modifications
		##					region 5 -eleRange 2121 2132 rayleigh 0.0 0.0 $a1 0.0;		# assign stiffness proportional damping to frame beams & columns w/out n modifications
		##					region 6 -eleRange 500000 599999 rayleigh 0.0 0.0 $a1 0.0;	# assign stiffness proportional damping to panel zone elements
		## 					rayleigh $a0 0.0 0.0 0.0;              						# assign mass proportional damping to structure (only assigns to nodes with mass)
		## updated 9 May 2012:	use "-rayleigh" instead of "rayleigh" when defining damping with the "region" command
		## 						use "region" command when defining mass proportional damping so that the stiffness proportional damping isn't canceled
		region 4 -eleRange 111 213 -rayleigh 0.0 0.0 $a1_mod 0.0;		# assign stiffness proportional damping to frame beams & columns w/ n modifications
		region 5 -eleRange 2121 2132 -rayleigh 0.0 0.0 $a1 0.0;			# assign stiffness proportional damping to frame beams & columns w/out n modifications
		region 6 -eleRange 500000 599999 -rayleigh 0.0 0.0 $a1 0.0;		# assign stiffness proportional damping to panel zone elements
		region 7 -node 1205 1305 2205 2305 -rayleigh $a0 0.0 0.0 0.0;	# assign mass proportional damping to structure (only assigns to nodes with mass)


	# define ground motion parameters
		set patternID 1;				# load pattern ID
		set GMdirection 1;				# ground motion direction (1 = x)
		set GMfile "NR94cnp.txt";		# ground motion filename
		set dt 0.01;					# timestep of input GM file
		set Scalefact 1.0;				# ground motion scaling factor
		set TotalNumberOfSteps 2495;	# number of steps in ground motion
		set GMtime [expr $dt*$TotalNumberOfSteps + 10.0];	# total time of ground motion + 10 sec of free vibration

	# define the acceleration series for the ground motion
		# syntax:  "Series -dt $timestep_of_record -filePath $filename_with_acc_history -factor $scale_record_by_this_amount
		set accelSeries "Series -dt $dt -filePath $GMfile -factor [expr $Scalefact*$g]";

	# create load pattern: apply acceleration to all fixed nodes with UniformExcitation
		# command: pattern UniformExcitation $patternID $GMdir -accel $timeSeriesID
		pattern UniformExcitation $patternID $GMdirection -accel $accelSeries;

	# define dynamic analysis parameters
		set dt_analysis 0.001;			# timestep of analysis
		wipeAnalysis;					# destroy all components of the Analysis object, i.e. any objects created with system, numberer, constraints, integrator, algorithm, and analysis commands
		constraints Plain;				# how it handles boundary conditions
		numberer RCM;					# renumber dof's to minimize band-width (optimization)
		system UmfPack;					# how to store and solve the system of equations in the analysis
		test NormDispIncr 1.0e-8 10;	# type of convergence criteria with tolerance, max iterations
		algorithm Newton;				# use Newton's solution algorithm: updates tangent stiffness at every iteration
		integrator Newmark 0.5 0.25;	# uses Newmark's average acceleration method to compute the time history
		analysis Transient;				# type of analysis: transient or static
		set NumSteps [expr round(($GMtime + 0.0)/$dt_analysis)];	# number of steps in analysis

	# perform the dynamic analysis and display whether analysis was successful
                puts "Starting analysis over $NumSteps steps"
                # progress create [expr $NumSteps/10]
                # foreach i [range 0 $NumSteps 10] {
		#     set ok [analyze 10 $dt_analysis];
                #     progress update
                # }
		set ok [analyze $NumSteps $dt_analysis];	# ok = 0 if analysis was completed
		if {$ok == 0} {
			puts "Dynamic analysis complete";
		} else {
			puts "Dynamic analysis did not converge";
		}

	# output time at end of analysis
		set currentTime [getTime];	# get current analysis time	(after dynamic analysis)
		puts "The current time is: $currentTime";
		wipe all;
}

wipe all;
