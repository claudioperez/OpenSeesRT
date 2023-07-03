# --------------------------------------------------------------------------------------------------
# build a section
#			Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source LibUnits.tcl;			# define units

# MATERIAL parameters -------------------------------------------------------------------
set IDconcU 1; 			# material ID tag -- unconfined cover concrete
set IDreinf 2; 				# material ID tag -- reinforcement
# nominal concrete compressive strength
set fc 		[expr -4.0*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec 		[expr 57*$ksi*sqrt(-$fc/$psi)];	# Concrete Elastic Modulus
# unconfined concrete
set fc1U 		$fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U	-0.003;			# strain at maximum strength of unconfined concrete
set fc2U 		[expr 0.2*$fc1U];		# ultimate stress
set eps2U	-0.01;			# strain at ultimate stress
set lambda 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
set ftU [expr -0.14*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness
# -----------
set Fy 		[expr 66.8*$ksi];		# STEEL yield stress
set Es		[expr 29000.*$ksi];		# modulus of steel
set Bs		0.005;			# strain-hardening ratio 
set R0 18;				# control the transition from elastic to plastic branches
set cR1 0.925;				# control the transition from elastic to plastic branches
set cR2 0.15;				# control the transition from elastic to plastic branches
uniaxialMaterial Concrete02 $IDconcU $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# build cover concrete (unconfined)
uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;				# build reinforcement material

# section GEOMETRY -------------------------------------------------------------
set HSec [expr 5.*$ft]; 		# Column Depth
set BSec [expr 3.*$ft];		# Column Width
set coverSec [expr 5.*$in];		# Column cover to reinforcing steel NA.
set numBarsSec 4;			# number of longitudinal-reinforcement bars in steel layer. (symmetric top & bot)
set barAreaSec [expr 1*$in2];	# area of longitudinal-reinforcement bars
set SecTag 1;			# set tag for symmetric section

# FIBER SECTION properties -------------------------------------------------------------
# symmetric section
#                        y
#                        ^
#                        |     
#             ---------------------     --   --
#             |   o     o     o    |     |    -- cover
#             |                       |     |
#             |                       |     |
#    z <--- |          +           |     H
#             |                       |     |
#             |                       |     |
#             |   o     o     o    |     |    -- cover
#             ---------------------     --   --
#             |-------- B --------|
#
# RC section: 
   set coverY [expr $HSec/2.0];	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coverZ [expr $BSec/2.0];	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coreY [expr $coverY-$coverSec ];	# The distance from the section z-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set coreZ [expr $coverZ-$coverSec ];	# The distance from the section y-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set nfY 16;			# number of fibers for concrete in y-direction
   set nfZ 4;			# number of fibers for concrete in z-direction
   section fiberSec $SecTag   {;	# Define the fiber section
	patch quadr $IDconcU $nfZ $nfY -$coverY $coverZ -$coverY -$coverZ $coverY -$coverZ $coverY $coverZ; 	# Define the concrete patch
	layer straight $IDreinf $numBarsSec $barAreaSec  $coreY $coreZ  $coreY -$coreZ;	# top layer reinforcement
	layer straight $IDreinf $numBarsSec $barAreaSec -$coreY $coreZ -$coreY -$coreZ;	# bottom layer reinfocement
    };	# end of fibersection definition

