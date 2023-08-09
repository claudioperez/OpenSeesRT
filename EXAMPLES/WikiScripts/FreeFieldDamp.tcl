###########################################################
#                                                         #
# Site response analysis of a layered soil column using   # 
# total stress analysis.  The finite rigidity of the      #
# underlying medium is considered through the use of a    # 
# viscous damper at the base.                             #
#                                                         #
#   Created by:  Chris McGann                             #
#                HyungSuk Shin                            #
#                Pedro Arduino                            #
#                Peter Mackenzie-Helnwein                 #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN and m  (unless specified)       #
#                                                         #
###########################################################

wipe

#-------------------------------------------------------------------------------------------
#  1. DEFINE GEOMETRY OF SOIL PROFILE
#-------------------------------------------------------------------------------------------

# number of soil layers
set totalLayer      3

# layer thicknesses (m)
set layerThick(3)   2
set layerThick(2)   8
set layerThick(1)   40

#-------------------------------------------------------------------------------------------
#  2. DEFINE MESH GEOMETRY
#-------------------------------------------------------------------------------------------

# define horizontal size of elements in all layers (m)
set sizeEleX        1.0

# define number of elements in each layer (in vertical direction)
set numEleY(3)      2
set numEleY(2)      8 
set numEleY(1)      40

# number of nodes in each layer (in vertical direction)
set numNodeY(3)     [expr 2*$numEleY(3)]
set numNodeY(2)     [expr 2*$numEleY(2)]
set numNodeY(1)     [expr 2*($numEleY(1) + 1)]


# define vertical size of elements in each layer (m)
for {set i 1} {$i <= $totalLayer} {incr i} {
    set sizeEleY($i) [expr $layerThick($i)/$numEleY($i)]
}

# total number of elements in vertical direction
set totalNumEleY    0
for {set i 1} {$i <= $totalLayer} {incr i 1} {
    set totalNumEleY [expr $totalNumEleY + $numEleY($i)]
}

# total number of nodes
set totalNumNode    [expr 2*($totalNumEleY + 1)] 

#-------------------------------------------------------------------------------------------
#  3. DEFINE NODES FOR SOIL ELEMENTS
#-------------------------------------------------------------------------------------------

# soil nodes are created in 2 dimensions, with 3 dof (2 translational, 1 porePressure)
model BasicBuilder -ndm 2 -ndf 2

set yCoord     0.0
set count      0
set nextLayer  0
# loop over layers
for {set i 1} {$i <= $totalLayer} {incr i 1} {
  # loop over nodes
    for {set j 1} {$j <= $numNodeY($i)} {incr j 2} {
        node    [expr $j+$nextLayer]     0.0         [expr $yCoord + $count*$sizeEleY($i)]
        node    [expr $j+$nextLayer+1]   $sizeEleX   [expr $yCoord + $count*$sizeEleY($i)]
       #puts "node    [expr $j+$nextLayer]     0.0         [expr $yCoord + $count*$sizeEleY($i)]"
       #puts "node    [expr $j+$nextLayer+1]   $sizeEleX   [expr $yCoord + $count*$sizeEleY($i)]"

        set count [expr $count+1]
    }
    set nextLayer [expr $nextLayer + $j-1]
}
puts "Finished creating all soil nodes..."

#-------------------------------------------------------------------------------------------
#  4. DEFINE DASHPOT NODES
#-------------------------------------------------------------------------------------------

node 2000 0.0 0.0
node 2001 0.0 0.0

puts "Finished creating dashpot nodes..."

#-------------------------------------------------------------------------------------------
#  5. DEFINE BOUNDARY CONDITIONS AND EQUAL DOF
#-------------------------------------------------------------------------------------------

# define fixity of base nodes
fix 1 0 1
fix 2 0 1

# define fixity of dashpot nodes
fix  2000  1 1
fix  2001  0 1

# define equal DOF for simple shear deformation of soil elements
for {set k 3} {$k <= $totalNumNode} {incr k 2} {
    equalDOF  $k  [expr $k+1]  1 2
}

# define equal DOF for dashpot and base soil nodes
equalDOF 1    2 1 
equalDOF 1 2001 1

puts "Finished creating all boundary conditions and equalDOF..."

#-------------------------------------------------------------------------------------------
#  6. DEFINE MATERIAL PROPERTIES OF SOIL LAYERS
#-------------------------------------------------------------------------------------------

# mass density of each layer (Mg/m^3)
set massDensity(3)  1.7
set massDensity(2)  1.7
set massDensity(1)  1.9

# shear wave velocity, Vs, for each layer (m/s)
set Vs(3)           450.0
set Vs(2)           180.0
set Vs(1)           510.0

# maximum shear modulus, Gmax, for each layer (kPa)
for {set k 1} {$k <= $totalLayer} {incr k 1} {
    set G($k)    [expr $massDensity($k)*$Vs($k)*$Vs($k)]
}

# poisson's ratio for soil
set nu              0.35

# reference bulk modulus (kPa)
for {set k 1} {$k <= $totalLayer} {incr k 1} {
    set bulk($k)    [expr (2*$G($k)*(1 + $nu))/(3*(1 - 2*$nu))]
}

# friction angle at peak shear stress for each layer (degrees)
set phi(3)          39
set phi(2)          30
set phi(1)          40

# peak shear strain for all layers
set gammaPeak       0.1

# reference pressure for definition of material parameters (kPa)
set refPress        80

# pressure dependency coefficient
set pressCoeff      0.5

# phase transformation angle for each layer (degrees)
set ptAngle(3)      27
set ptAngle(2)      29
set ptAngle(1)      27

# contraction coefficient for each layer
set contract(3)     0.05
set contract(2)     0.21
set contract(1)     0.05

# first rate of dilation coefficient for each layer 
set dilate1(3)      0.6
set dilate1(2)      0.0
set dilate1(1)      0.6

# second rate of dilation coefficient for each layer 
set dilate2(3)      3.0
set dilate2(2)      0.0
set dilate2(1)      3.0

# liquefaction parameters (no water in this example --> no liquefaction potential)
set liq1            0.0
set liq2            0.0
set liq3            0.0

# inital void ratio for each layer
set void(3)         0.55
set void(2)         0.85
set void(1)         0.55

# number of yield surfaces
set numSurf         20

# critical state line parameters
set critState1      0.9
set critState2      0.02
set critState3      0.7
set pAtm            101

#-------------------------------------------------------------------------------------------
#  7. DEFINE SOIL MATERIALS
#-------------------------------------------------------------------------------------------

for {set i 1} {$i <= $totalLayer} {incr i 1} {
    nDMaterial PressureDependMultiYield $i 2 $massDensity($i) $G($i) $bulk($i) $phi($i) $gammaPeak \
                                        $refPress $pressCoeff $ptAngle($i) $contract($i) $dilate1($i) \
                                        $dilate2($i) $liq1 $liq2 $liq3 $numSurf $void($i) $critState1 \
                                        $critState2 $critState3 $pAtm 
}
puts "Finished creating all soil materials..."

#-------------------------------------------------------------------------------------------
#  7. DEFINE SOIL ELEMENTS
#-------------------------------------------------------------------------------------------

set wgtX 0.0

for {set i 1} {$i <= $totalLayer} {incr i 1} {
    set wgtY($i)  [expr -9.81*$massDensity($i)]
}

set eleTag  1
# loop over layers
for {set i 1} {$i <= $totalLayer} {incr i 1} {
  # loop over elements
    for {set j 1} {$j <= $numEleY($i)} {incr j 1} {

        set nI  [expr 2*$eleTag - 1] 
        set nJ  [expr $nI + 1]
        set nK  [expr $nI + 3]
        set nL  [expr $nI + 2]

        element quad $eleTag $nI $nJ $nK $nL 1.0 "PlaneStrain" $i 0.0 0.0 $wgtX $wgtY($i)
        #puts "quad $eleTag $nI $nJ $nK $nL 1.0 PlaneStrain $i 0.0 0.0 $wgtX $wgtY($i)"

        set eleTag  [expr $eleTag + 1]
    }
}
updateMaterialStage -material 1 -stage 0
updateMaterialStage -material 2 -stage 0
updateMaterialStage -material 3 -stage 0

puts "Finished creating all soil elements..."

#-------------------------------------------------------------------------------------------
#  8. DEFINE MATERIAL AND ELEMENTS FOR VISCOUS DAMPERS
#-------------------------------------------------------------------------------------------

# bedrock shear wave velocity (m/s)
set bedVS  914.4
# bedrock mass density (Mg/m^3)
set bedDen 2.4

# dashpot coefficient 
set mC     [expr $bedDen*$bedVS]

# material
uniaxialMaterial Viscous 4000 $mC 1

# elements
element zeroLength 5000 2000 2001 -mat 4000 -dir 1

puts "Finished creating dashpot material and element..."

#-------------------------------------------------------------------------------------------
#  8. CREATE RECORDERS
#-------------------------------------------------------------------------------------------

# record nodal displacments, velocities, and accelerations at each time step
recorder Node -file displacement.out -time -nodeRange 1 $totalNumNode -dof 1 2  disp
recorder Node -file velocity.out     -time -nodeRange 1 $totalNumNode -dof 1 2  vel
recorder Node -file acceleration.out -time -nodeRange 1 $totalNumNode -dof 1 2  accel

# record stress and strain at each gauss point in the soil elements
recorder Element -file stress1.out   -time -eleRange  1   $totalNumEleY   material 1 stress
recorder Element -file stress2.out   -time -eleRange  1   $totalNumEleY   material 2 stress
recorder Element -file stress3.out   -time -eleRange  1   $totalNumEleY   material 3 stress
recorder Element -file stress4.out   -time -eleRange  1   $totalNumEleY   material 4 stress

recorder Element -file strain1.out   -time -eleRange  1   $totalNumEleY   material 1 strain
recorder Element -file strain2.out   -time -eleRange  1   $totalNumEleY   material 2 strain
recorder Element -file strain3.out   -time -eleRange  1   $totalNumEleY   material 3 strain
recorder Element -file strain4.out   -time -eleRange  1   $totalNumEleY   material 4 strain

# real time display recorder for visualization during analysis
# # recorder display "OpenSees Real Time" 10 10 700 700 -wipe
# prp            0 0 100
# vup            0 1 0
# vpn            0 0 1
# display        1 3 100

puts "Finished creating all recorders..."

#-------------------------------------------------------------------------------------------
#  9. APPLY GRAVITY LOADING
#-------------------------------------------------------------------------------------------

# Newmark parameters
set gamma   0.6
set beta    [expr pow($gamma+0.5, 2)/4]

# analysis objects
constraints Transformation
test        NormDispIncr 1e-5 30 1
algorithm   ModifiedNewton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
analysis    Transient

analyze     10 5.0e2

puts "Finished with elastic gravity analysis..."

#-------------------------------------------------------------------------------------------
#  10. UPDATE SOIL MATERIAL STAGE
#-------------------------------------------------------------------------------------------

# change to plastic behavior
for {set i 1} {$i <= $totalLayer} {incr i} {
    updateMaterialStage -material $i -stage 1
}

analyze    40 5.0e2

puts "Finished with plastic gravity analysis..."

#-------------------------------------------------------------------------------------------
#  11. HORIZONTAL DYNAMIC ANALYSIS  
#-------------------------------------------------------------------------------------------

setTime 0.0
wipeAnalysis

set gFactor 1.0

# information from force history file
set forceFile forceHistory.out
# time step increment
set dT 0.02
# number of data points
set nSteps 4184

# timeseries object for force history
set mSeries "Path -dt $dT -filePath $forceFile -factor $gFactor"

# loading object
pattern Plain 10 $mSeries {
load 1 1.0 0.0 
}
puts "Dynamic loading created..."

# analysis objects
constraints Penalty 1.0e16 1.0e16
test        NormDispIncr 1e-3 30 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
analysis    Transient

analyze     $nSteps  $dT

puts "Finished with dynamic analysis..."

