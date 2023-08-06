#-*-opensees-*-
#YamamotoBiaxialHDR
#2-directional loading for EESD 2012;41:1815-1860,eqe2161

#begin
wipe
set outDir ./Output/
file mkdir $outDir

#---------------------------------------- model
model basic -ndm 3 -ndf 6 ;#BasicBuilder

#---------------------------------------- node,fix
node 1  0.0 0.0 0.0
node 2  0.0 0.0 1.0

fix  1   1 1 1 1 1 1
fix  2   0 0 1 1 1 1

#---------------------------------------- material, element
uniaxialMaterial Elastic 99  1.0e9

#                             iN jN Tp   DDo   DDi     Hr
element YamamotoBiaxialHDR  1  1  2  1  1.30  0.03  0.261  -orient 0 0 1 1 0 0 -coRS 1.0 1.0; #local(x,y,z)=Global(Z,X,Y)
element twoNodeLink         2  1  2  -mat 99 99  -dir 2 3  -orient 0 0 1 1 0 0 ;#control local-y&z deformation

#---------------------------------------- recorder
recorder Element -file ./$outDir/YamamotoBiaxialHDR_def.txt -time -ele 1 localDisplacement
recorder Element -file ./$outDir/YamamotoBiaxialHDR_frc.txt -time -ele 1 localForce

#---------------------------------------- load
timeSeries Path 1 -dt 1.0 -filePath ./YamamotoBiaxialHDR_input_X.tcl
timeSeries Path 2 -dt 1.0 -filePath ./YamamotoBiaxialHDR_input_Y.tcl

pattern Plain 1 1 {
    load 2  [expr 0.001*1.0e9] 0 0 0 0 0
}
pattern Plain 2 2 {
    load 2  0 [expr 0.001*1.0e9] 0 0 0 0;#bi-directional
#    load 2  0 [expr 0.000*1.0e9] 0 0 0 0;#uni-directional
}

#---------------------------------------- analysis
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 6 ;# tol,numIter
algorithm Newton

#---------------------------------------- run
#displacement control with stiff spring
integrator LoadControl 1.0;
analysis Static
analyze 901

#---------------------------------------- print
print ele 1

#end
remove recorders

