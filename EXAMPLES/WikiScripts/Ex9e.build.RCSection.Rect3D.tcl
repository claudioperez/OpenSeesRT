# --------------------------------------------------------------------------------------------------
# build a section
#			Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source LibUnits.tcl;			# define units

# MATERIAL parameters -------------------------------------------------------------------
set IDconcCore 1; 				# material ID tag -- confined core concrete
set IDconcCover 2; 				# material ID tag -- unconfined cover concrete
set IDreinf 3; 				# material ID tag -- reinforcement
# nominal concrete compressive strength
set fc 		[expr -4.0*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec 		[expr 57*$ksi*sqrt(-$fc/$psi)];	# Concrete Elastic Modulus
# confined concrete
set Kfc 		1.3;			# ratio of confined to unconfined concrete strength
set fc1C 		[expr $Kfc*$fc];		# CONFINED concrete (mander model), maximum stress
set eps1C	[expr 2.*$fc1C/$Ec];	# strain at maximum stress 
set fc2C 		[expr 0.2*$fc1C];		# ultimate stress
set eps2C 	[expr 5*$eps1C];		# strain at ultimate stress 
# unconfined concrete
set fc1U 		$fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U	-0.003;			# strain at maximum strength of unconfined concrete
set fc2U 		[expr 0.2*$fc1U];		# ultimate stress
set eps2U	-0.01;			# strain at ultimate stress
set lambda 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
set ftC [expr -0.14*$fc1C];		# tensile strength +tension
set ftU [expr -0.14*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness
# -----------
set Fy 		[expr 66.8*$ksi];		# STEEL yield stress
set Es		[expr 29000.*$ksi];		# modulus of steel
set Bs		0.01;			# strain-hardening ratio 
set R0 18;				# control the transition from elastic to plastic branches
set cR1 0.925;				# control the transition from elastic to plastic branches
set cR2 0.15;				# control the transition from elastic to plastic branches
uniaxialMaterial Concrete02 $IDconcCore $fc1C $eps1C $fc2C $eps2C $lambda $ftC $Ets;	# build core concrete (confined)
uniaxialMaterial Concrete02 $IDconcCover $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# build cover concrete (unconfined)
uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;				# build reinforcement material

# section GEOMETRY -------------------------------------------------------------
set HSec [expr 5.*$ft]; 		# Column Depth
set BSec [expr 3.*$ft];		# Column Width
set coverH [expr 5.*$in];		# Column cover to reinforcing steel NA, parallel to H
set coverB [expr 3.5*$in];		# Column cover to reinforcing steel NA, parallel to B
set numBarsTop 16;		# number of longitudinal-reinforcement bars in steel layer. -- top
set numBarsBot 16;		# number of longitudinal-reinforcement bars in steel layer. -- bot
set numBarsIntTot 6;			# number of longitudinal-reinforcement bars in steel layer. -- total intermediate skin reinforcement, symm about y-axis
set barAreaTop [expr 2.25*$in2];	# area of longitudinal-reinforcement bars -- top
set barAreaBot [expr 2.25*$in2];	# area of longitudinal-reinforcement bars -- bot
set barAreaInt [expr 2.25*$in2];	# area of longitudinal-reinforcement bars -- intermediate skin reinf
set SecTag 1;			# set tag for symmetric section

# FIBER SECTION properties -------------------------------------------------------------
#
#                        y
#                        ^
#                        |     
#             ---------------------    --   --
#             |   o     o     o    |     |    -- coverH
#             |                      |     |
#             |   o            o    |     |
#    z <--- |          +          |     Hsec
#             |   o            o    |     |
#             |                      |     |
#             |   o o o o o o    |     |    -- coverH
#             ---------------------    --   --
#             |-------Bsec------|
#             |---| coverB  |---|
#
#                       y
#                       ^
#                       |    
#             ---------------------
#             |\      cover        /|
#             | \------Top------/ |
#             |c|                   |c|
#             |o|                   |o|
#  z <-----|v|       core      |v|  Hsec
#             |e|                   |e|
#             |r|                    |r|
#             | /-------Bot------\ |
#             |/      cover        \|
#             ---------------------
#                       Bsec
#
# Notes
#    The core concrete ends at the NA of the reinforcement
#    The center of the section is at (0,0) in the local axis system

set coverY [expr $HSec/2.0];	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
set coverZ [expr $BSec/2.0];	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
set coreY [expr $coverY-$coverH];	# The distance from the section z-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
set coreZ [expr $coverZ-$coverB];	# The distance from the section y-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concreteset nfY 16;			# number of fibers for concrete in y-direction
set nfY 16;			# number of fibers for concrete in y-direction
set nfZ 4;				# number of fibers for concrete in z-direction
set numBarsInt [expr $numBarsIntTot/2];	# number of intermediate bars per side
section fiberSec $SecTag     {;	# Define the fiber section
	patch quadr $IDconcCore $nfZ $nfY -$coreY $coreZ -$coreY -$coreZ $coreY -$coreZ $coreY $coreZ; 	# Define the core patch
	patch quadr $IDconcCover 1 $nfY -$coverY $coverZ -$coreY $coreZ $coreY $coreZ $coverY $coverZ;	# Define the four cover patches
	patch quadr $IDconcCover 1 $nfY -$coreY -$coreZ -$coverY -$coverZ $coverY -$coverZ $coreY -$coreZ
	patch quadr $IDconcCover $nfZ 1 -$coverY $coverZ -$coverY -$coverZ -$coreY -$coreZ -$coreY $coreZ
	patch quadr $IDconcCover $nfZ 1 $coreY $coreZ $coreY -$coreZ $coverY -$coverZ $coverY $coverZ
	layer straight $IDreinf $numBarsInt $barAreaInt  -$coreY $coreZ $coreY $coreZ;	# intermediate skin reinf. +z
	layer straight $IDreinf $numBarsInt $barAreaInt  -$coreY -$coreZ $coreY -$coreZ;	# intermediate skin reinf. -z
	layer straight $IDreinf $numBarsTop $barAreaTop $coreY $coreZ $coreY -$coreZ;	# top layer reinfocement
	layer straight $IDreinf $numBarsBot $barAreaBot  -$coreY $coreZ  -$coreY -$coreZ;	# bottom layer reinforcement
};	# end of fibersection definition


# assign torsional Stiffness for 3D Model
set SecTagTorsion 99;		# ID tag for torsional section behavior
set SecTag3D 3;			# ID tag for combined behavior for 3D model
uniaxialMaterial Elastic $SecTagTorsion $Ubig;	# define elastic torsional stiffness
section Aggregator $SecTag3D $SecTagTorsion T -section $SecTag;	# combine section properties
