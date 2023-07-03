# CenterColShearSpring.tcl
# Defines Shear Spring 
# Units: kip, in
# KJE, Feb 2003


# slopes of shear spring backbone
set rigidSlope 1700; #Values when using zero-length spring (G*Av/L)
set negSlope -8;

# residual shear capacity
set Vr 3.0;

# strengths for initial response 
set Vi1 25.0;
set Vi2 30.0;
set Vi3 45.0;

# stiffness of unloading slope for flexural component
set kf 24.7 ;# measured off hysteresis plot from analysis

# define limit surface using shear drift model 
#                     tag        eleTag 
# rho   f'c   b   h    d   Fsw
# Kdeg Fres defType forType nodeI nodeJ dof perpDirn 
limitCurve Shear $shearCurveTag $bcTag\
    0.0018 3517.0 9.0 9.0 7.75 11.87\
    $kf    $Vr   2       0        1     4    1    2  0.0

# define HystereticMaterial
uniaxialMaterial LimitState $shearTag\
    $Vi1 [expr $Vi1/$rigidSlope] $Vi2 [expr $Vi2/$rigidSlope] $Vi3 [expr $Vi3/$rigidSlope]\
    [expr -$Vi1] [expr -$Vi1/$rigidSlope]  [expr -$Vi2] [expr -$Vi2/$rigidSlope] [expr -$Vi3] [expr -$Vi3/$rigidSlope]\
    $pinchX $pinchY $damage1 $damage2 $beta $shearCurveTag 2 0



