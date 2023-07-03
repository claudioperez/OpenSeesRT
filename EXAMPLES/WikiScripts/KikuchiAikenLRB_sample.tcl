#-*-opensees-*-
#KikuchiAikenLRB sample

#----- begin
wipe

#---------------------------------------- model, node, fix
model basic -ndm 1 -ndf 1

node 1 0.0
node 2 0.0

fix 1 1
fix 2 0

#---------------------------------------- material, element
#type=1, Ar=0.7540[m^2], Hr=0.200[m], Gr=0.4e6[N/m^2], Ap=0.0314[m^2], Tp=8.33e6[N/m^2], Alph=0.588e6[N/m^2], Beta=13.0
uniaxialMaterial KikuchiAikenLRB 1  1 0.7540 0.200 0.4e6 0.0314 8.33e6 0.588e6 13.0

element zeroLength 1 1 2 -mat 1 -dir 1

#---------------------------------------- recorder
recorder Element -file KikuchiAikenLRB_output.txt -ele 1 deformationsANDforces

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

integrator DisplacementControl 2 1 0.002
analyze 100
print node 2

integrator DisplacementControl 2 1 -0.002
analyze 200

integrator DisplacementControl 2 1 0.002
analyze 300
print node 2

integrator DisplacementControl 2 1 -0.002
analyze 400

integrator DisplacementControl 2 1 0.002
analyze 500
print node 2

integrator DisplacementControl 2 1 -0.002
analyze 600

integrator DisplacementControl 2 1 0.002
analyze 700
print node 2

integrator DisplacementControl 2 1 -0.002
analyze 800

integrator DisplacementControl 2 1 0.002
analyze 400

#----- end
remove recorders
