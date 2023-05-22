# This example implements a slight modification of the verification test from reference #1.
#

model BasicBuilder -ndm 3 -ndf 3

set x 30. ;     #Set to example from paper x = 30, 60, 80, 100. Will not work for x=0.01, system ill-conditioned.

node 1 0.0 0.0 90.0
node 2 [expr $x/2] 0.0 40.0
node 3 $x 60 30.

fix 1 1 1 1 
fix 2 0 1 0 
fix 3 1 1 1

set w3  -0.00001
set E  3.e7 
set A  1.
set L0  100. 
set alfa  6.5e-6 
set cambiodetemp  100.
set rho [expr $w3 / 9.81]

set errorTol 1e-6
set NSubSteps 20

element CatenaryCable 1 1 2 $w3 $E $A [expr $L0/2] $alfa $cambiodetemp $rho $errorTol $NSubSteps  0
element CatenaryCable 2 2 3 $w3 $E $A [expr $L0/2] $alfa $cambiodetemp $rho $errorTol $NSubSteps 0

set NSteps 10
timeSeries Linear 1 -factor 1

pattern Plain 2 1 {
    eleLoad -ele 1 2 -type -beamUniform 0. 0. -1
}


recorder Node -file out/disp.txt -time -nodeRange 1 3 -dof 1 2 3 disp
recorder Element -file out/forces.txt -time -eleRange 1 2 force

system FullGeneral
constraints Plain
numberer Plain
test NormDispIncr 1.0e-5 100 1
integrator LoadControl [expr 1.0/$NSteps]
algorithm Newton
analysis Static

analyze $NSteps

print -node 2

