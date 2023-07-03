#-*-opensees-*-
#AxialSp sample

#----- begin
wipe

#---------------------------------------- model, node, fix
model basic -ndm 1 -ndf 1

node 1 0.0
node 2 0.0

fix 1 1
fix 2 0

#---------------------------------------- material, element
#compressive modulus = 1000e6[Pa]
#yield stress = 1e6[Pa] (tension), -25e6[Pa] (compression)
#reduction rate = 0.50 (tensile), 0.01 (tensile yielding), 0.50(compressive yeilding)
#target point Stress = -5e6[Pa]
uniaxialMaterial AxialSp 1 1000e6 1e6 -25e6 0.50 0.01 0.50 -5e6

element zeroLength 1 1 2 -mat 1 -dir 1

#---------------------------------------- recorder
recorder Element -file AxialSp_output.txt -ele 1 deformationsANDforces

#---------------------------------------- load
timeSeries Linear 1

pattern Plain 1 1 {
    load 2 1.0
}

#---------------------------------------- analysis
constraints Transformation
numberer RCM
system BandGeneral
test NormUnbalance 1.0e-5 10
algorithm Newton

#---------------------------------------- run
integrator DisplacementControl 2 1 0.0
analysis Static
analyze 1

integrator DisplacementControl 2 1 0.001
analyze 10
print node 2

integrator DisplacementControl 2 1 -0.001
analyze 20

integrator DisplacementControl 2 1 0.001
analyze 30
print node 2

integrator DisplacementControl 2 1 -0.001
analyze 40
print node 2

integrator DisplacementControl 2 1 0.001
analyze 50
print node 2

integrator DisplacementControl 2 1 -0.001
analyze 60
print node 2

integrator DisplacementControl 2 1 0.001
analyze 30
print node 2

#----- end
remove recorders
