# --------------------------------------------------------------------------------------------------
# build a section
#		Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source LibUnits.tcl;			# define units

# MATERIAL parameters -------------------------------------------------------------------
# define MATERIAL properties ----------------------------------------
set Fy [expr 60.0*$ksi]
set Es [expr 29000*$ksi];		# Steel Young's Modulus
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
set Hiso 0
set Hkin 1000
set matIDhard 1
uniaxialMaterial Hardening  $matIDhard $Es $Fy   $Hiso  $Hkin

# Structural-Steel W-section properties -------------------------------------------------------------------
set SecTag 1
set WSec W27x114

# from Steel Manuals:
# in × lb/ft 	Area (in2) 	d (in) 	bf (in) 	tf (in) 	tw (in) 	Ixx (in4) 	Iyy (in4)
# W27x114  	33.5 		27.29 	10.07 	0.93 	0.57 	4090 	159
set d [expr 27.29*$in];	# nominal depth
set tw [expr 0.57*$in];	# web thickness
set bf [expr 10.07*$in];	# flange width
set tf [expr 0.93*$in];	# flange thickness
set nfdw 16;	# number of fibers along web depth 
set nftw 4;	# number of fibers along web thickness
set nfbf 16;	# number of fibers along flange width (you want this many in a bi-directional loading)
set nftf 4;	# number of fibers along flange thickness
  
  set dw [expr $d - 2 * $tf]
  set y1 [expr -$d/2]
  set y2 [expr -$dw/2]
  set y3 [expr  $dw/2]
  set y4 [expr  $d/2]
  
  set z1 [expr -$bf/2]
  set z2 [expr -$tw/2]
  set z3 [expr  $tw/2]
  set z4 [expr  $bf/2]
  
  #                           
  section fiberSec  $SecTag {
     #                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
     patch quadr  $matIDhard $nfbf $nftf   $y1 $z4   $y1 $z1   $y2 $z1   $y2 $z4
     patch quadr  $matIDhard $nftw $nfdw   $y2 $z3   $y2 $z2   $y3 $z2   $y3 $z3
     patch quadr  $matIDhard $nfbf $nftf   $y3 $z4   $y3 $z1   $y4 $z1   $y4 $z4
  }

# assign torsional Stiffness for 3D Model
set SecTagTorsion 99;		# ID tag for torsional section behavior
set SecTag3D 3;			# ID tag for combined behavior for 3D model
uniaxialMaterial Elastic $SecTagTorsion $Ubig;	# define elastic torsional stiffness
section Aggregator $SecTag3D $SecTagTorsion T -section $SecTag;	# combine section properties
