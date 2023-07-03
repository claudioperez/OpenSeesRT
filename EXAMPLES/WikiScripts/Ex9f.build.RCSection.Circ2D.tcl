# --------------------------------------------------------------------------------------------------
# build a section
#		Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
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
set ftC [expr -0.14*$fc1C];			# tensile strength +tension
set ftU [expr -0.14*$fc1U];			# tensile strength +tension
set Ets [expr $ftU/0.002];			# tension softening stiffness
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
set DSec [expr 5.*$ft]; 		# Column Diameter
set coverSec [expr 5.*$in];		# Column cover to reinforcing steel NA.
set numBarsSec 16;		# number of uniformly-distributed longitudinal-reinforcement bars
set barAreaSec [expr 2.25*$in2];	# area of longitudinal-reinforcement bars
set SecTag 1;			# set tag for symmetric section

# Generate a circular reinforced concrete section
# with one layer of steel evenly distributed around the perimeter and a confined core.
# confined core.
#		by:  Michael H. Scott, 2003
# 
#
# Notes
#    The center of the reinforcing bars are placed at the inner radius
#    The core concrete ends at the inner radius (same as reinforcing bars)
#    The reinforcing bars are all the same size
#    The center of the section is at (0,0) in the local axis system
#    Zero degrees is along section y-axis
# 
set ri 0.0;			# inner radius of the section, only for hollow sections
set ro [expr $DSec/2];	# overall (outer) radius of the section
set nfCoreR 8;		# number of radial divisions in the core (number of "rings")
set nfCoreT 8;		# number of theta divisions in the core (number of "wedges")
set nfCoverR 4;		# number of radial divisions in the cover
set nfCoverT 8;		# number of theta divisions in the cover

# Define the fiber section
section fiberSec $SecTag  {
	set rc [expr $ro-$coverSec];					# Core radius
	patch circ $IDconcCore $nfCoreT $nfCoreR 0 0 $ri $rc 0 360;		# Define the core patch
	patch circ $IDconcCover $nfCoverT $nfCoverR 0 0 $rc $ro 0 360;	# Define the cover patch
	set theta [expr 360.0/$numBarsSec];		# Determine angle increment between bars
	layer circ $IDreinf $numBarsSec $barAreaSec 0 0 $rc $theta 360;	# Define the reinforcing layer
}