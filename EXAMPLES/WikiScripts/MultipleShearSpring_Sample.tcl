#-*-opensees-*-
#MultipleShearSpring sample

#----- begin
wipe

#---------------------------------------- model, node, fix
model basic -ndm 3 -ndf 6

node 1   0.0 0.0 0.0
node 2   0.0 0.0 1.0

fix 1   1 1 1 1 1 1
fix 2   0 0 1 1 1 1

#---------------------------------------- material, element
uniaxialMaterial Steel01 1   25 100 0.05
uniaxialMaterial Elastic 99   1.0e9

element multipleShearSpring 1   1 2   16 -mat 1 -lim 0.5 -orient 0 0 1 1 0 0 ;#nSpring=16, local(x,y,z)=Global(Z,X,Y)
element twoNodeLink 99   1 2   -mat 99 99  -dir 2 3 -orient 0 0 1 1 0 0 ;#control local-y&z deformation of MSS

#---------------------------------------- recorder
recorder Element -file ./MultipleShearSpring_output_deformation.txt -time -ele 1 basicDeformation
recorder Element -file ./MultipleShearSpring_output_force.txt -time -ele 1 basicForce

#---------------------------------------- load
timeSeries Path 1 -dt 1.0 -filePath ./MultipleShearSpring_input_X.tcl
timeSeries Path 2 -dt 1.0 -filePath ./MultipleShearSpring_input_Y.tcl

pattern Plain 1 1 {
    load 2   1.0e9 0 0 0 0 0 ;
}
pattern Plain 2 2 {
    load 2   0 1.0e9 0 0 0 0 ;
}

#---------------------------------------- analysis
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 6
algorithm Newton

#---------------------------------------- run
#displacement control with stiff spring
integrator LoadControl 1.0
analysis Static
analyze 1

integrator LoadControl 1.0
analyze 250
print node 2

integrator LoadControl 1.0
analyze 250
print node 2

integrator LoadControl 1.0
analyze 250
print node 2

#----- end
remove recorders
