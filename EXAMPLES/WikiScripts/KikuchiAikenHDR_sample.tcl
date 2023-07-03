#-*-opensees-*-
#KikuchiAikenHDR sample

#----- begin
wipe

#---------------------------------------- model, node, fix
model basic -ndm 1 -ndf 1

node 1 0.0
node 2 0.0

fix 1 1
fix 2 0

#---------------------------------------- material, element
uniaxialMaterial KikuchiAikenHDR 1 X0.6 1.0 1.0 ;#type=X0.6, Ar=1.0[m^2], Hr=1.0[m]

element zeroLength 1 1 2 -mat 1 -dir 1

#---------------------------------------- recorder
recorder Element -file KikuchiAikenHDR_output.txt -ele 1 deformationsANDforces

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

integrator DisplacementControl 2 1 0.01
analyze 100
print node 2

integrator DisplacementControl 2 1 -0.01
analyze 200

integrator DisplacementControl 2 1 0.01
analyze 300
print node 2

integrator DisplacementControl 2 1 -0.01
analyze 400

integrator DisplacementControl 2 1 0.01
analyze 500
print node 2

integrator DisplacementControl 2 1 -0.01
analyze 600

integrator DisplacementControl 2 1 0.01
analyze 700
print node 2

integrator DisplacementControl 2 1 -0.01
analyze 800

integrator DisplacementControl 2 1 0.01
analyze 400

#----- end
remove recorders
