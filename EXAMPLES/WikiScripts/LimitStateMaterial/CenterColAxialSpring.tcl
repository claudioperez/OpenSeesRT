#CenterColAxialSpring.tcl

# Defines Axial Spring 
# Units: kip, in
# KJE, Feb 2003

set Fsw [expr 11.87/1.0]

# Axial elastic stiffness (only one axial spring used)
set axialElasticSlope [expr 99.0*$Ec*$A/$L];  #99 times more rigid than column
set axialNegSlope -90.0;  #slope = axial load/axial displ}

# residual capacity 
set Pr 5.0;

# axial loads used for setting initial elastic slope 
set P1 65.0;
set P2 75.0;
set P3 85.0;


# define limit surface
#                     tag       eleTag   Fsw  Kdeg          Fres defType forType nodeI nodeJ dof perpDir delta eleRemove
limitCurve Axial $axialCurveTag $bcTag  $Fsw $axialNegSlope $Pr   2       2      1     4     1   2       0.0   0

# define LimitStateMaterial 
uniaxialMaterial LimitState $axialFailTag\
    $P1 [expr $P1/$axialElasticSlope]  $P2 [expr $P2/$axialElasticSlope] $P3 [expr $P3/$axialElasticSlope]\
    [expr -$P1] [expr -$P1/$axialElasticSlope] [expr -$P2] [expr -$P2/$axialElasticSlope] [expr -$P3] [expr -$P3/$axialElasticSlope]\
    0.5 0.5 0.0 0.0 0.0 $axialCurveTag 1

