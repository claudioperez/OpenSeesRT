###########################################################
#                                                         #
# 1D consolidation using fully coupled, u-p formulation.  #
#  Cohesive soil. Option of single or double drainage.    #
#                                                         #
#   Created by:  HyungSuk Shin                            #
#                Chris McGann                             #
#                Pedro Arduino                            #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN, m, and sec (unless specified)  #
#                                                         #
###########################################################

wipe

# specify desired drainage type; 1 --> single drainage, 2 --> double drainage
set drainage      1

#-------------------------------------------------------------------------------------------
#  1. DEFINE MESH GEOMETRY
#-------------------------------------------------------------------------------------------

# number of elements in horizontal direction
set numEleX       1
# size of elements in horizontal direction (m)
set sizeEleX      1.0
# number of nodes in horizontal direction
set numNodeX      [expr 2*$numEleX + 1]

# number of elements in vertical direction
set numEleY       10
# size of elements in vertical direction (m)
set sizeEleY      1.0
# number of nodes in vertical direction 
set numNodeY      [expr 2*$numEleY + 1]

# total number of nodes
set totalNumNode  [expr 3*$numNodeY]
# total number of elements
set totalNumEle   [expr $numEleY*$numEleX]

#-------------------------------------------------------------------------------------------
#  2. DEFINE CORNER NODES
#-------------------------------------------------------------------------------------------

# corner nodes are created in 2D with 3 DOF
model BasicBuilder -ndm 2 -ndf 3

# loop in horizontal direction
for {set i 1} {$i <= $numNodeX} {incr i 2} {
    # loop in vertical direction
    for {set j 1} {$j <= $numNodeY} {incr j 2} {

        set xCoord  [expr ($i-1)*$sizeEleX/2] 
        set yCoord  [expr ($j-1)*$sizeEleX/2] 
        set nodeNum [expr $i + ($j-1)*$numNodeX] 
        
        node  $nodeNum  $xCoord  $yCoord
        puts "node  $nodeNum  $xCoord  $yCoord"
  }
}
puts "Finished creating corner soil nodes..."

#-------------------------------------------------------------------------------------------
#  3. DEFINE BOUNDARY CONDITIONS FOR CORNER NODES
#-------------------------------------------------------------------------------------------

# boundary conditions for base nodes depending upon selected drainage type
switch $drainage {
    1 {
        fix  1  1 1 0
        fix  3  1 1 0
    }
    2 {
        fix  1  1 1 1
        fix  3  1 1 1
    }
}


# fix horizontal displacement and porepressure DOF for nodes at the top
fix [expr $totalNumNode-2]  1 0 1
fix $totalNumNode           1 0 1

# fix horizontal displacement for all remaining nodes
for {set i 1} {$i <= [expr 3*$numNodeY-6]} {incr i 6} {
    if {$i == 1 || $i == 3} continue
    puts "fix  $i  1 0 0"
    fix  $i  1 0 0
    puts "fix  [expr $i+2]  1 0 0"
    fix  [expr $i+2]  1 0 0
}

puts "Finished creating boundary conditions for corner nodes..."

#-------------------------------------------------------------------------------------------
#  4. DEFINE INTERIOR NODES
#-------------------------------------------------------------------------------------------

# interior nodes are created in 2D with 2 DOF
model BasicBuilder -ndm 2 -ndf 2

# central column of nodes
set xCoord  [expr $sizeEleX/2]

for {set j 1} {$j <= $numNodeY} {incr j 1} {
   
    set yCoord  [expr ($j-1)*$sizeEleX/2]
    set nodeNum [expr 3*$j - 1] 

    node  $nodeNum  $xCoord  $yCoord 
    puts "node  $nodeNum  $xCoord  $yCoord"
}

# interior nodes on the element edges
for {set j 1} {$j <= $numEleY} {incr j 1} {

    set yCoord   [expr $j - 0.5]
    set nodeNumL [expr 6*$j - 2]
    set nodeNumR [expr $nodeNumL + 2]
    
    node  $nodeNumL  0.0  $yCoord
    node  $nodeNumR  $sizeEleX  $yCoord
    puts "node  $nodeNumL  0.0  $yCoord"
    puts "node  $nodeNumR  $sizeEleX  $yCoord"
}
puts "Finished creating all soil nodes..."

#-------------------------------------------------------------------------------------------
#  5. DEFINE BOUNDARY CONDITIONS FOR INTERIOR NODES
#-------------------------------------------------------------------------------------------

# fix displacement DOF for node at the base
fix  2  1 1

# fix horizontal displacement of nodes on edge of elements
for {set i 1} {$i <= $numEleY} {incr i 1} {

    set nodeNumL [expr 6*$i - 2]
    set nodeNumR [expr $nodeNumL + 2]

    fix  $nodeNumL  1 0
    fix  $nodeNumR  1 0
    puts "fix  $nodeNumL  1 0"
    puts "fix  $nodeNumR  1 0"
}

# define equalDOF for surface nodes
puts "equalDOF [expr $totalNumNode-1] [expr $totalNumNode-2]  1 2"
puts "equalDOF [expr $totalNumNode-1] [expr $totalNumNode]    1 2"
equalDOF [expr $totalNumNode-1] [expr $totalNumNode-2]  1 2
equalDOF [expr $totalNumNode-1] [expr $totalNumNode]    1 2

puts "Finished creating all boundary conditions and equalDOF..."

#-------------------------------------------------------------------------------------------
#  6. DEFINE MATERIAL PROPERTIES OF SOIL
#-------------------------------------------------------------------------------------------

# saturated mass density (Mg/m3)
set satDensity  2.3
# fluid mass density (Mg/m3)
set H2ODensity  1.0

# shear modulus (kPa)
set shear       2.5e4
# bulk modulus (kPa)
set bulk        6.2e5
# undrained bulk modulus (kPa)
set uBulk       2.2e5

# vertical permeability, input value in m/s, units are modified to suit constitutive model
set mPermV      5.0e-5
set vPerm       [expr $mPermV/9.81/$H2ODensity]

# horizontal permeability, input value in m/s, units are modified to suit constitutive model
set mPermH      5.0e-5
set hPerm       [expr $mPermH/9.81/$H2ODensity]

# cohesion (kPa)
set cohesion    45.0

# friction angle (degrees)
set phi         0.0

# peak shear strain 
set gammaPeak   0.1

# reference pressure (kPa)
set refPress    80.0

# pressure dependency coefficient
set pressCoeff  0.0

# number of yield surfaces
set numSurf     22

#-------------------------------------------------------------------------------------------
#  7. DEFINE SOIL MATERIALS
#-------------------------------------------------------------------------------------------

set matTag 1
nDMaterial PressureIndependMultiYield $matTag 2 $satDensity $shear $bulk $cohesion $gammaPeak \
                                      $phi $refPress $pressCoeff $numSurf

puts "Finished creating all soil materials..."

#-------------------------------------------------------------------------------------------
#  8. DEFINE SOIL ELEMENTS
#-------------------------------------------------------------------------------------------

set thick 1.0
set wgtX  0.0
set wgtY  -9.81

for {set j 1} {$j <= $numEleY} {incr j 1} {
    
    set nI  [expr 6*$j - 5]
    set nJ  [expr $nI + 2]
    set nK  [expr $nI + 8]
    set nL  [expr $nI + 6]
    set nM  [expr $nI + 1]
    set nN  [expr $nI + 5]
    set nP  [expr $nI + 7]
    set nQ  [expr $nI + 3]
    set nR  [expr $nI + 4]

    element 9_4_QuadUP $j $nI $nJ $nK $nL $nM $nN $nP $nQ $nR \
                       $thick $matTag $bulk $H2ODensity $hPerm $vPerm $wgtX $wgtY 
} 

# update material stage to zero for elastic behavior
updateMaterialStage -material $matTag -stage 0

puts "Finished creating all soil elements..."

#-------------------------------------------------------------------------------------------
#  8. APPLY GRAVITY LOADING
#-------------------------------------------------------------------------------------------

# Newmark parameters
set gamma   0.5
set beta    0.25

# analysis objects
constraints Transformation
test        NormDispIncr 1e-5 30 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
analysis    Transient

analyze     10 5.0e2

puts "Finished with elastic gravity analysis..."

#-------------------------------------------------------------------------------------------
#  9. UPDATE SOIL MATERIAL STAGE
#-------------------------------------------------------------------------------------------

# update material stage to consider elastoplastic behavior
updateMaterialStage -material $matTag -stage 1

analyze    40 5.0e2

puts "Finished with plastic gravity analysis..."


#-------------------------------------------------------------------------------------------
#  10. CREATE RECORDERS
#-------------------------------------------------------------------------------------------

# record nodal displacments, porewater pressures, and accelerations 
# recorder Node -file Output/displacement.out -time -node all -dof 1 2  disp
# recorder Node -file Output/porepressure.out -time -node all -dof 3    vel

# record stress and strain at third Gauss point in the elements
recorder Element -file out/stress.out  -time -eleRange  1   $numEleY  material 9 stress
recorder Element -file out/strain.out  -time -eleRange  1   $numEleY  material 9 strain

puts "Finished creating all recorders..."

#-------------------------------------------------------------------------------------------
#  11. CONSOLIDATION ANALYSIS  
#-------------------------------------------------------------------------------------------

setTime 0.0
wipeAnalysis

set overburden -1000.0

# constant loading applied at the surface
pattern Plain 1 "Constant" {
    load  [expr $totalNumNode-1]  0.0 $overburden
}
puts "Consolidation loading created..."

# analysis objects
constraints Penalty 1.0e16 1.0e16
test        NormDispIncr 1e-4 30 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
rayleigh    0.5 0.2 0.0 0.0
analysis    Transient

analyze     800 0.25

puts "Finished with consolidation analysis..."

wipe  

