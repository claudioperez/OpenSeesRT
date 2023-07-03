# --------------------------------------------------------------------------------------------------
# Moment-Curvature analysis of section
#		Silvia Mazzoni & Frank McKenna, 2006
#

# define procedure
source MomentCurvature3D.tcl

# set AXIAL LOAD --------------------------------------------------------
set P [expr -1800*$kip];	# + Tension, - Compression

# set maximum Curvature:
set Ku [expr 0.01/$in];
set numIncr 100;	# Number of analysis increments to maximum curvature (default=100)
# Call the section analysis procedure
MomentCurvature3D $SecTag3D $P $Ku $numIncr