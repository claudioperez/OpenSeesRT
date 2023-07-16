###########################################################
#                                                         #
# 3D site response analysis of a soil deposit on an       #
# elastic half-space.  Shaking is applied in a single     #
# plane, and the site has a slope out of the plane of     #
# shaking.  The finite rigidity of the underlying medium  #
# is considered through the use of a viscous damper at    #
# the base of the soil column.                            #
#                                                         #
#   Created by:  Chris McGann                             #
#                Pedro Arduino                            #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN and m   (unless specified)      #
#                                                         #
###########################################################

wipe
set outDir ./out/
file mkdir $outDir
#-----------------------------------------------------------------------------------------
#  1. DEFINE ANALYSIS PARAMETERS
#-----------------------------------------------------------------------------------------

#---SOIL GEOMETRY
# thickness of soil profile (m)
set soilThick       30.0
# grade of slope (%)
set grade           2.0
# number of soil layers
set numLayers       3
# layer thicknessess (m)
set layerThick(3)   2.0
set layerThick(2)   8.0
set layerThick(1)   20.0
# depth of water table
set waterTable     2.0

# define layer boundaries
set layerBound(1) $layerThick(1)
for {set i 2} {$i <= $numLayers} {incr i 1} {
    set layerBound($i) [expr $layerBound([expr $i-1])+$layerThick($i)]
}

#---MESH GEOMETRY
# number of elements in horizontal direction
set nElemX  1
set nElemZ  1
# number of nodes in horizontal direction
set nNodeX  [expr 2*$nElemX+1]
# horizontal element size (m)
set sElemX  2.0
set sElemZ  2.0

# number of elements in vertical direction for each layer
set nElemY(3)  4
set nElemY(2)  16
set nElemY(1)  40
# total number of elements in vertical direction
set nElemT     60
# vertical element size in each layer
for {set i 1} {$i <=$numLayers} {incr i 1} {
    set sElemY($i) [expr $layerThick($i)/$nElemY($i)]
    puts "size:  $sElemY($i)"
}

# number of nodes in vertical direction in each layer
set nNodeT 0
for {set k 1} {$k < $numLayers} {incr k 1} {
    set nNodeY($k)  [expr 4*$nElemY($k)]
    puts "number of nodes in layer $k: $nNodeY($k)"
    set nNodeT  [expr $nNodeT + $nNodeY($k)]
}
set nNodeY($numLayers) [expr 4*($nElemY($numLayers) + 1)]
puts "number of nodes in layer $numLayers: $nNodeY($numLayers)"
set nNodeT  [expr $nNodeT + $nNodeY($numLayers)]
puts "total number of nodes: $nNodeT"

#-----------------------------------------------------------------------------------------
#  2. CREATE SOIL NODES AND BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 3 -ndf 4

set yCoord  0.0 
set count   0
set gwt     1
set waterHeight [expr $soilThick-$waterTable]
set nodesInfo [open $outDir/nodesInfo.dat w]
# loop over layers
for {set k 1} {$k <= $numLayers} {incr k 1} {
	# loop over nodes
	for {set j 1} {$j <= $nNodeY($k)} {incr j 4} {
		node  [expr $j+$count]    0.0      $yCoord  0.0
		node  [expr $j+$count+1]  0.0      $yCoord  $sElemZ
		node  [expr $j+$count+2]  $sElemX  $yCoord  $sElemZ
		node  [expr $j+$count+3]  $sElemX  $yCoord  0.0

		puts $nodesInfo "[expr $j+$count]    0.0      $yCoord  0.0"
		puts $nodesInfo "[expr $j+$count+1]  0.0      $yCoord  $sElemZ"
		puts $nodesInfo "[expr $j+$count+2]  $sElemX  $yCoord  $sElemZ"
		puts $nodesInfo "[expr $j+$count+3]  $sElemX  $yCoord  0.0"

		# designate nodes above water table
        if {$yCoord>=$waterHeight} {
            set dryNode($gwt) [expr $j+$count]
            set gwt [expr $gwt+1]
			set dryNode($gwt) [expr $j+$count+1]
            set gwt [expr $gwt+1]
			set dryNode($gwt) [expr $j+$count+2]
            set gwt [expr $gwt+1]
			set dryNode($gwt) [expr $j+$count+3]
            set gwt [expr $gwt+1]
        }

		set yCoord  [expr $yCoord + $sElemY($k)]
	}
	set count  [expr $count + $nNodeY($k)]
}
close $nodesInfo
puts "Finished creating all soil nodes..."

# define boundary conditions for nodes at base of column
fix 1  0 1 1 0
fix 2  0 1 1 0
fix 3  0 1 1 0
fix 4  0 1 1 0
equalDOF  1 2 1
equalDOF  1 3 1
equalDOF  1 4 1

# define periodic boundary conditions for remaining nodes
set count  0
for {set k 5} {$k <= [expr $nNodeT]} {incr k 4} {
	equalDOF  $k  [expr $k+1]  1 2 3
	equalDOF  $k  [expr $k+2]  1 2 3
	equalDOF  $k  [expr $k+3]  1 2 3
	puts "equalDOF  $k  [expr $k+1]  1 2 3"
	puts "equalDOF  $k  [expr $k+2]  1 2 3"
	puts "equalDOF  $k  [expr $k+3]  1 2 3"
}

# define pore pressure boundaries for nodes above water table
for {set i 1} {$i < $gwt} {incr i 1} {
    fix $dryNode($i)  0 0 0 1
}
puts "Finished creating all soil boundary conditions..."

#-----------------------------------------------------------------------------------------
#  4. CREATE SOIL MATERIALS
#-----------------------------------------------------------------------------------------

set slope [expr atan($grade/100.0)]
set g     -9.81

nDMaterial PressureDependMultiYield02 3 3 1.8 9.0e4 2.2e5 32 0.1 \
                                      101.0 0.5 26 0.067 0.23 0.06 \
                                      0.27 20 5.0 3.0 1.0 \
                                      0.0 0.77 0.9 0.02 0.7 101.0
set xWgt(3)  0.0
set yWgt(3)  [expr $g*cos($slope)]
set zWgt(3)  [expr $g*sin($slope)]
set uBulk(3) 2.2e6
set xPerm(3) 1.0e-2
set yPerm(3) 1.0e-2
set zPerm(3) 1.0e-2
set eVoid(3) 0.77

nDMaterial PressureDependMultiYield02 2 3 2.24 9.0e4 2.2e5 32 0.1 \
                                      101.0 0.5 26 0.067 0.23 0.06 \
                                      0.27 20 5.0 3.0 1.0 \
                                      0.0 0.77 0.9 0.02 0.7 101.0
set xWgt(2)  0.0
set yWgt(2)  [expr $g*cos($slope)]
set zWgt(2)  [expr $g*sin($slope)]
set uBulk(2) 2.2e6
set xPerm(2) 1.0e-5
set yPerm(2) 1.0e-5
set zPerm(2) 1.0e-5
set eVoid(2) 0.77

nDMaterial PressureDependMultiYield02 1 3 2.45 1.3e5 2.6e5 39 0.1 \
                                      101.0 0.5 26 0.010 0.0 0.35 \
                                      0.0 20 5.0 3.0 1.0 \
                                      0.0 0.47 0.9 0.02 0.7 101.0
set xWgt(1)  0.0
set yWgt(1)  [expr $g*cos($slope)]
set zWgt(1)  [expr $g*sin($slope)]
set uBulk(1) 2.2e6
set xPerm(1) 1.0e-2
set yPerm(1) 1.0e-2
set zPerm(1) 1.0e-2
set eVoid(1) 0.47

puts "Finished creating all soil materials..."

#-----------------------------------------------------------------------------------------
#  5. CREATE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------

set alpha 1.5e-6
set count 0
set elemInfo [open $outDir/elementInfo.dat w]
# loop over layers 
for {set k 1} {$k <= $numLayers} {incr k 1} {
    # loop over elements
    for {set j 1} {$j <= $nElemY($k)} {incr j 1} {

        set nI  [expr 4*($j+$count) - 3] 
        set nJ  [expr $nI + 1]
        set nK  [expr $nI + 2]
        set nL  [expr $nI + 3]
		set nM  [expr $nI + 4]
		set nN  [expr $nI + 5]
		set nO  [expr $nI + 6]
		set nP  [expr $nI + 7]

		element SSPbrickUP [expr $j+$count] $nI $nJ $nK $nL $nM $nN $nO $nP $k $uBulk($k) 1.0 1.0 1.0 1.0 $eVoid($k) $alpha $xWgt($k) $yWgt($k) $zWgt($k)
		puts $elemInfo "[expr $j+$count] $nI $nJ $nK $nL $nM $nN $nO $nP"
    }
    set count [expr $count + $nElemY($k)]
}
close $elemInfo
puts "Finished creating all soil elements..."

#-----------------------------------------------------------------------------------------
#  6. LYSMER DASHPOTS
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 3 -ndf 3

# define dashpot nodes
set dashF [expr $nNodeT + 1]
set dashX [expr $nNodeT + 2]
set dashZ [expr $nNodeT + 3]

node $dashF  0.0 0.0 0.0
node $dashX  0.0 0.0 0.0
node $dashZ  0.0 0.0 0.0

# define fixities for dashpot nodes
fix $dashF  1 1 1
fix $dashX  0 1 1
fix $dashZ  1 1 0

# link dashpots with soil column
equalDOF  1 $dashX 1
equalDOF  1 $dashZ 3

# dashpot material
set colArea       [expr $sElemX*$sElemZ]
set rockVS        700.0
set rockDen       2.5
set dashpotCoeff  [expr $rockVS*$rockDen*$colArea]
uniaxialMaterial Viscous [expr $numLayers+1] $dashpotCoeff 1

# dashpot elements
element zeroLength [expr $nElemT+1]  $dashF $dashX -mat [expr $numLayers+1] -dir 1
element zeroLength [expr $nElemT+2]  $dashF $dashZ -mat [expr $numLayers+1] -dir 3

puts "Finished creating Lysmer dashpots..."

#-----------------------------------------------------------------------------------------
#  8. CREATE GID FLAVIA.MSH FILE FOR POSTPROCESSING
#-----------------------------------------------------------------------------------------
 
set meshFile [open freeField3D.flavia.msh w]
puts $meshFile "MESH ffBrick dimension 3 ElemType Hexahedra Nnode 8"
puts $meshFile "Coordinates"
puts $meshFile "#node_number   coord_x   coord_y   coord_z"
set yCoord  0.0 
set count   0
# loop over layers
for {set k 1} {$k <= $numLayers} {incr k 1} {
	# loop over nodes
	for {set j 1} {$j <= $nNodeY($k)} {incr j 4} {
		puts $meshFile  "[expr $j+$count]    0.0      $yCoord  0.0"
		puts $meshFile  "[expr $j+$count+1]  0.0      $yCoord  $sElemZ"
		puts $meshFile  "[expr $j+$count+2]  $sElemX  $yCoord  $sElemZ"
		puts $meshFile  "[expr $j+$count+3]  $sElemX  $yCoord  0.0"

		set yCoord  [expr $yCoord + $sElemY($k)]
	}
	set count  [expr $count + $nNodeY($k)]
}
puts $meshFile "end coordinates"
puts $meshFile "Elements"
puts $meshFile "# element   node1   node2   node3   node4   node5   node6   node7   node8"
set count 0
# loop over layers 
for {set k 1} {$k <= $numLayers} {incr k 1} {
    # loop over elements
    for {set j 1} {$j <= $nElemY($k)} {incr j 1} {

        set nI  [expr 4*($j+$count) - 3] 
        set nJ  [expr $nI + 1]
        set nK  [expr $nI + 2]
        set nL  [expr $nI + 3]
		set nM  [expr $nI + 4]
		set nN  [expr $nI + 5]
		set nO  [expr $nI + 6]
		set nP  [expr $nI + 7]

        puts $meshFile  "[expr $j+$count] $nI $nJ $nK $nL $nM $nN $nO $nP"
    }
    set count [expr $count + $nElemY($k)]
}
puts $meshFile "end elements"
close $meshFile

#-----------------------------------------------------------------------------------------
#  7. GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------

## record nodal displacments, velocities, and accelerations at each time step
recorder Node -file $outDir/Gdisplacement.out -time  -nodeRange 1 $nNodeT -dof 1 2 3 disp
recorder Node -file $outDir/Gvelocity.out     -time  -nodeRange 1 $nNodeT -dof 1 2 3 vel
recorder Node -file $outDir/Gacceleration.out -time  -nodeRange 1 $nNodeT -dof 1 2 3 accel
recorder Node -file $outDir/GporePressure.out -time  -nodeRange 1 $nNodeT -dof 4 vel

# record stress and strain at each gauss point in the soil elements
recorder Element -file $outDir/Gstress.out   -time  -eleRange  1   $nElemT   stress
recorder Element -file $outDir/Gstrain.out   -time  -eleRange  1   $nElemT  strain

puts "Finished creating gravity recorders..."

#-----------------------------------------------------------------------------------------
#  9. DEFINE ANALYSIS PARAMETERS
#-----------------------------------------------------------------------------------------

#---GROUND MOTION PARAMETERS
# time step in ground motion record
set motionDT     0.005
# number of steps in ground motion record
set motionSteps  7990

#---RAYLEIGH DAMPING PARAMETERS
set pi      3.141592654
# damping ratio
set damp    0.02
# lower frequency
set omega1  [expr 2*$pi*0.2]
# upper frequency
set omega2  [expr 2*$pi*20]
# damping coefficients
set a0      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set a1      [expr 2*$damp/($omega1 + $omega2)]
puts "damping coefficients: a_0 = $a0;  a_1 = $a1"

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
# maximum shear wave velocity (m/s)
set vsMax       250.0
# duration of ground motion (s)
set duration    [expr $motionDT*$motionSteps]
# minimum element size
set minSize $sElemY(1)
for {set i 2} {$i <= $numLayers} {incr i 1} {
    if {$sElemY($i) < $minSize} {
        set minSize $sElemY($i)
    }
}
# trial analysis time step
set kTrial      [expr $minSize/(pow($vsMax,0.5))]
# define time step and number of steps for analysis
if { $motionDT <= $kTrial } {
    set nSteps  $motionSteps
    set dT      $motionDT
} else {
    set nSteps  [expr int(floor($duration/$kTrial)+1)]
    set dT      [expr $duration/$nSteps]
}
puts "number of steps in analysis: $nSteps"
puts "analysis time step: $dT"

#---ANALYSIS PARAMETERS
# Newmark parameters
set gamma  0.5
set beta   0.25

#-----------------------------------------------------------------------------------------
#  7. GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 3 -ndf 4

# update materials to consider elastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 0
}

constraints Penalty 1.e14 1.e14
test        NormDispIncr 1e-5 30 1
algorithm   Newton
numberer    Plain
system      SparseGeneral
integrator  Newmark $gamma $beta 
analysis    Transient

set startT  [clock seconds]
analyze     20 5.0e2

puts "Finished with elastic gravity analysis..."

# update materials to consider plastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 1
}

# plastic gravity loading
analyze     40 5.0e-2

puts "Finished with plastic gravity analysis..."

#-----------------------------------------------------------------------------------------
#  11. UPDATE ELEMENT PERMEABILITY VALUES FOR POST-GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

# layer 3
setParameter -value $xPerm(3) -eleRange 57 60 xPerm
setParameter -value $yPerm(3) -eleRange 57 60 yPerm
setParameter -value $zPerm(3) -eleRange 57 60 zPerm
# layer 2
setParameter -value $xPerm(2) -eleRange 41 56 xPerm
setParameter -value $yPerm(2) -eleRange 41 56 yPerm
setParameter -value $zPerm(2) -eleRange 41 56 zPerm
# layer 1
setParameter -value $xPerm(1) -eleRange 1 40 xPerm
setParameter -value $yPerm(1) -eleRange 1 40 yPerm
setParameter -value $zPerm(1) -eleRange 1 40 zPerm

puts "Finished updating permeabilities for dynamic analysis..."

#-----------------------------------------------------------------------------------------
#  8. CREATE POST-GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------
 
# reset time and analysis
setTime 0.0
wipeAnalysis
remove recorders
 
# recorder time step
set recDT  [expr 10*$motionDT]

## record nodal displacments, velocities, and accelerations at each time step
recorder Node -file displacement.out -time -dT $recDT  -nodeRange 1 $nNodeT -dof 1 2 3 disp
recorder Node -file velocity.out     -time -dT $recDT  -nodeRange 1 $nNodeT -dof 1 2 3 vel
recorder Node -file acceleration.out -time -dT $recDT  -nodeRange 1 $nNodeT -dof 1 2 3 accel
recorder Node -file porePressure.out -time -dT $recDT  -nodeRange 1 $nNodeT -dof 4 vel

# record stress and strain at each gauss point in the soil elements
recorder Element -file stress.out   -time -dT $recDT -eleRange  1   $nElemT   stress
recorder Element -file strain.out   -time -dT $recDT -eleRange  1   $nElemT  strain

puts "Finished creating post-gravity recorders..."

#-----------------------------------------------------------------------------------------
#  9. DYNAMIC ANALYSIS
#-----------------------------------------------------------------------------------------

# define velocity time history file
set velocityFile yerbaNSvelocity.out
 
# timeseries object for force history
set mSeries "Path -dt $motionDT -filePath $velocityFile -factor $dashpotCoeff"
 
# loading object
pattern Plain 10 $mSeries {
    load 1  1.0 0.0 0.0 0.0
}
puts "Dynamic loading created..."
 
constraints Penalty 1.e14 1.e14
test        NormDispIncr 1.0e-3 55 1
algorithm   KrylovNewton
numberer    Plain
system      SparseGeneral
integrator  Newmark $gamma $beta
#rayleigh    $a0 $a1 0.0 0.0
analysis    Transient

# perform analysis with timestep reduction loop
set ok [analyze $nSteps  $dT]
 
# if analysis fails, reduce timestep and continue with analysis
if {$ok != 0} {
    puts "did not converge, reducing time step"
    set curTime  [getTime]
    set mTime $curTime
    puts "curTime: $curTime"
    set curStep  [expr $curTime/$dT]
    puts "curStep: $curStep"
    set rStep  [expr ($nSteps-$curStep)*2.0]
    set remStep  [expr int(($nSteps-$curStep)*2.0)]
    puts "remStep: $remStep"
    set dT       [expr $dT/2.0]
    puts "dT: $dT"
 
    set ok [analyze  $remStep  $dT]
 
    # if analysis fails again, reduce timestep and continue with analysis
    if {$ok != 0} {
        puts "did not converge, reducing time step"
        set curTime  [getTime]
        puts "curTime: $curTime"
        set curStep  [expr ($curTime-$mTime)/$dT]
        puts "curStep: $curStep"
        set remStep  [expr int(($rStep-$curStep)*2.0)]
        puts "remStep: $remStep"
        set dT       [expr $dT/2.0]
        puts "dT: $dT"
 
        analyze  $remStep  $dT
    }
}
set endT    [clock seconds]
puts "Finished with dynamic analysis..."
puts "Analysis execution time: [expr $endT-$startT] seconds"

wipe
