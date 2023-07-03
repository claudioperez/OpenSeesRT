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
set SecTagFlex 2;			# assign a tag number to the column flexural behavior
set SecTagAxial 3;			# assign a tag number to the column axial behavior	
set SecTag 1;			# assign a tag number to the column section tag

# COLUMN section
# calculated stiffness parameters
set EASec $Ubig;				# assign large value to axial stiffness
set MySec [expr 130000*$kip*$in];		# yield moment
set PhiYSec [expr 0.65e-4/$in];		# yield curvature
set EICrack [expr $MySec/$PhiYSec];		# cracked section inertia
set b 0.01 ;				# strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
uniaxialMaterial Steel01 $SecTagFlex $MySec $EICrack $b; 		# bilinear behavior for flexural moment-curvature
uniaxialMaterial Elastic $SecTagAxial $EASec;			# this is not used as a material, this is an axial-force-strain response

# assign torsional Stiffness for 3D Model
set SecTagTorsion 99;		# ID tag for torsional section behavior
set SecTag3D 5;			# ID tag for combined behavior for 3D model
uniaxialMaterial Elastic $SecTagTorsion $Ubig;	# define elastic torsional stiffness
section Aggregator $SecTag3D $SecTagAxial P $SecTagFlex Mz $SecTagTorsion T;	# combine section properties
