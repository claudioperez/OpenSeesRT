#-*-opensees-*-
#KikuchiBearing sample

#begin
wipe

#---------------------------------------- model, node, fix
model basic -ndm 3 -ndf 6

node 1   0.0 0.0 0.0
node 2   0.0 0.0 0.43856
node 3   0.0 0.0 0.43856

fix 1   1 1 1 1 1 1
fix 2   0 1 0 1 1 1
fix 3   1 1 1 1 1 1

#---------------------------------------- material, element
#diameter (rubber) = 1.016 [m]
#diameter (lead plug)= 0.200 [m]
#thickness (rubber) = 0.008 * 40 [m]
#thickness (steel) = 0.00304 * 39 [m]
#shear modulus (rubber) = 0.49e6 [N/m^2]
#bulk modulus (rubber) = 1960e6 [N/m^2]
#yield stress (lead plug) = 8.33e6 [N/m^2]
#shear modulus (lead plug) = 0.588e6 [N/m^2]
#ratio of initial stiffness to yielding stiffness = 13.0
#totalRubber = 0.320 [m]
#totalHeight = 0.43856 [m]
#area (rubber) = 0.7793 [m^2]
#area (lead plug) = 0.0314 [m^2]
#initial compression modulus = 1013e6 [N/m^2]
#limDisp = totalRubber * 25% = 0.080 [m]
#lambda = 1.016/0.008*sqrt(3*0.49e6/1960e6) = 3.478

uniaxialMaterial KikuchiAikenLRB 1   1 0.7793 0.320 0.49e6 0.0314 8.33e6 0.588e6 13.0
uniaxialMaterial AxialSp 2   1013e6 1e6 -100e6 1.00 0.01 0.50 0e6
uniaxialMaterial Elastic 99   1e12

#orient: local(x,y,z)=Global(Z,X,Y)
#element KikuchiBearing 1   1 2  -shape round -size 1.016 0.320 -nMSS 8 -matMSS 1 -nMNS 30 -matMNS 2 -orient 0 0 1 1 0 0 -limDisp 0.080 -noPDInput -noTilt ;#<<case 1>>
#element KikuchiBearing 1   1 2  -shape round -size 1.016 0.320 -nMSS 8 -matMSS 1 -nMNS 30 -matMNS 2 -orient 0 0 1 1 0 0 -limDisp 0.080 ;#<<case 2>>
element KikuchiBearing 1   1 2  -shape round -size 1.016 0.320 -nMSS 8 -matMSS 1 -nMNS 30 -matMNS 2 -orient 0 0 1 1 0 0 -limDisp 0.080 -lambda 3.478 ;#<<case 3>>

element zeroLength 99   3 2  -mat 99 -dir 2 -orient 0 0 1 1 0 0 ;#control local-y deformation

#---------------------------------------- recorder
recorder Element -file KikuchiBearing_output_deformation.txt -time -ele 1 basicDeformation
recorder Element -file KikuchiBearing_output_force.txt -time -ele 1 basicForce

#---------------------------------------- load
timeSeries Path 1 -dt 1.0 -filePath KikuchiBearing_input_Z.tcl
timeSeries Path 2 -dt 1.0 -filePath KikuchiBearing_input_X.tcl

pattern Plain 1 1 {
    load 2  0 0 -14000e3 0 0 0 ;#axial stress 17.3MPa
}

pattern Plain 2 2 {
    load 2  0.320e12 0 0 0 0 0 ;#shear strain 100%, 200%
}

#---------------------------------------- analysis
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-4 50
algorithm ModifiedNewton

#---------------------------------------- run
#displacement control with stiff spring
integrator LoadControl 0.002
analysis Static
analyze 1000
print node 2

integrator LoadControl 0.002
analyze 500
print node 2

integrator LoadControl 0.002
analyze 1500

integrator LoadControl 0.002
analyze 500
print node 2

integrator LoadControl 0.002
analyze 1500

#----- end
remove recorders
