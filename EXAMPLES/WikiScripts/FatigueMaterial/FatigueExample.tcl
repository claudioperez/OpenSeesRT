# Run an example scrript to test the fatigue material model

## Define model and loads
model      BasicBuilder   -ndm 2       -ndf 2
node 1 0.0 0.0
node 2 0.0 0.0
fix 1  1 1 
fix 2  0 1 
uniaxialMaterial Steel01 1 60.0 29800.0 0.003
uniaxialMaterial Fatigue 2 1
element zeroLength 1 1 2 -mat 2 -dir 1
pattern Plain  1   "Linear" {
    #    nd       FX   
    load 2    1.0 0.0 0.0
}

## Recorders
recorder Element -file "Damage.out" -time -ele 1 material 1 damage
recorder Element -file "StressStrain.out" \
    -time -ele 1 material 1 stressANDstrain

## Set analysis parameters
test EnergyIncr 1.0e-8  200    0
algorithm Newton 
system UmfPack
numberer RCM
constraints Plain
analysis Static

## Source the displacement history, and initialize analysis parameters
source RandomStrainHstory.tcl
set LoopLength [array size disp]
set h 1
set controlNode 2
set currentDisp [nodeDisp $controlNode 1 ]
puts [format " \n STARTING  DISPLACEMENT = %5.3f \n" $currentDisp]

## Run the static cyclic analysis
while {$h < $LoopLength} {
    set controlNodeDisp [nodeDisp $controlNode 1 ]   
    set dU [expr  (1.8*($disp($h) - $controlNodeDisp))/100.0]
    integrator DisplacementControl $controlNode  1 $dU 1 $dU $dU 
    set ok [ analyze 100]
    set h [expr $h + 1 ]    
}
