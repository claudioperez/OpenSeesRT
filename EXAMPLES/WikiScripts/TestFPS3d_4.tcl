# File: TestFPS3d_4.tcl
#
# $Revision: $
# $Date: $
# $URL: $
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Created: 02/09
# Revision: A
#
# Purpose: this file tests the 3D flatSliderBearing or the
# singleFPBearing element. It models an isolated five story  
# one bay building and the bearing element has finite length.
# It also tests the different friction models.

# ------------------------------
# Start of model generation
# ------------------------------
# Create ModelBuilder (with three-dimensions and 6 DOF/node)
model BasicBuilder -ndm 3 -ndf 6

# Define geometry for model
# -------------------------
set g [expr 32.174*12.0]
set P 3.0
set mass [expr $P/$g]
#    tag   xCrd   yCrd   zCrd        mass
node   1    0.0    0.0    0.0
node   2    0.0    0.0   10.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   3    0.0    0.0  154.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   4    0.0    0.0  298.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   5    0.0    0.0  442.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   6    0.0    0.0  586.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   7    0.0    0.0  730.0  -mass $mass $mass $mass 0.0 0.0 0.0
node   8  144.0    0.0    0.0
node   9  144.0    0.0   10.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  10  144.0    0.0  154.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  11  144.0    0.0  298.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  12  144.0    0.0  442.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  13  144.0    0.0  586.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  14  144.0    0.0  730.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  15    0.0  144.0    0.0
node  16    0.0  144.0   10.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  17    0.0  144.0  154.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  18    0.0  144.0  298.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  19    0.0  144.0  442.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  20    0.0  144.0  586.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  21    0.0  144.0  730.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  22  144.0  144.0    0.0
node  23  144.0  144.0   10.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  24  144.0  144.0  154.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  25  144.0  144.0  298.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  26  144.0  144.0  442.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  27  144.0  144.0  586.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  28  144.0  144.0  730.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  29   72.0   72.0   10.0;   # master node for 1st floor diaphragm
node  30   72.0   72.0  154.0;   # master node for 2nd floor diaphragm
node  31   72.0   72.0  298.0;   # master node for 3rd floor diaphragm
node  32   72.0   72.0  442.0;   # master node for 4th floor diaphragm
node  33   72.0   72.0  586.0;   # master node for 5th floor diaphragm
node  34   72.0   72.0  730.0;   # master node for 6th floor diaphragm

# Set the boundary conditions
#   tag DX DY DZ RX RY RZ
fix  1  1  1  1  1  1  1
fix  8  1  1  1  1  1  1
fix 15  1  1  1  1  1  1
fix 22  1  1  1  1  1  1
fix 29  0  0  1  1  1  0
fix 30  0  0  1  1  1  0
fix 31  0  0  1  1  1  0
fix 32  0  0  1  1  1  0
fix 33  0  0  1  1  1  0
fix 34  0  0  1  1  1  0

# Set the multi-point constraints 
# rigidDiaphragm perpDir mNodeTag sNodeTags 
rigidDiaphragm      3      29     2  9 16 23 
rigidDiaphragm      3      30     3 10 17 24 
rigidDiaphragm      3      31     4 11 18 25 
rigidDiaphragm      3      32     5 12 19 26 
rigidDiaphragm      3      33     6 13 20 27 
rigidDiaphragm      3      34     7 14 21 28 

# Define material models
# ----------------------
set mv [expr 6.0*$mass] 
set kv 7500.0
set zetaVertical 0.02
set cv [expr 2.0*$zetaVertical*sqrt($kv*$mv)]
uniaxialMaterial Elastic 1 $kv $cv
uniaxialMaterial Elastic 2 0.0

# Define friction model for FP elements
# -------------------------------------
# frictionModel Coulomb tag mu
frictionModel Coulomb 1 0.163

# frictionModel VDependent tag muSlow muFast transRate
#frictionModel VDependent 1 0.085 0.163 0.77

# frictionModel VPDependent tag muSlow muFast0 A deltaMu alpha transRate
#frictionModel VPDependent 1 0.085 0.163 7.0686 0.05 0.08 0.77

# Define geometric transformations
# --------------------------------
# geomTransf type tag vec_xz
geomTransf Linear  1  1  0 0
geomTransf Linear  2  0  1 0
geomTransf Linear  3  0 -1 0

# Define elements
# ---------------
# element flatSliderBearing eleTag NodeI NodeJ frnMdlTag kInit -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element flatSliderBearing 1  1  2 1 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element flatSliderBearing 2  8  9 1 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element flatSliderBearing 3 15 16 1 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element flatSliderBearing 4 22 23 1 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0

# element singleFPBearing eleTag NodeI NodeJ frnMdlTag Reff kInit -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
element singleFPBearing 1  1  2 1 34.68 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0 -iter 100 1E-12
element singleFPBearing 2  8  9 1 34.68 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0 -iter 100 1E-12
element singleFPBearing 3 15 16 1 34.68 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0 -iter 100 1E-12
element singleFPBearing 4 22 23 1 34.68 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0 -iter 100 1E-12

# element RJWatsonEqsBearing eleTag NodeI NodeJ frnMdlTag kInit k2 k3 mu -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element RJWatsonEqsBearing 1  1  2 1 250.0 0.519 0.0 3.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element RJWatsonEqsBearing 2  8  9 1 250.0 0.519 0.0 3.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element RJWatsonEqsBearing 3 15 16 1 250.0 0.519 0.0 3.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0
#element RJWatsonEqsBearing 4 22 23 1 250.0 0.519 0.0 3.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0

# element elasticBeamColumn eleTag NodeI NodeJ A E G J Iy Iz geoTranTag <-mass massDens> 
element elasticBeamColumn  5  2  3 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn  6  3  4 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn  7  4  5 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn  8  5  6 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn  9  6  7 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 10  9 10 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 11 10 11 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 12 11 12 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 13 12 13 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 14 13 14 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 15 16 17 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 16 17 18 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 17 18 19 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 18 19 20 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 19 20 21 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 20 23 24 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 21 24 25 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 22 25 26 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 23 26 27 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 24 27 28 20.0 29000.0 11154.0 100.0 400.0 400.0 2
element elasticBeamColumn 25  2  9 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 26 16 23 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 27  2 16 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 28  9 23 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 29  3 10 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 30 17 24 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 31  3 17 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 32 10 24 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 33  4 11 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 34 18 25 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 35  4 18 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 36 11 25 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 37  5 12 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 38 19 26 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 39  5 19 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 40 12 26 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 41  6 13 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 42 20 27 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 43  6 20 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 44 13 27 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 45  7 14 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 46 21 28 20.0 29000.0 11154.0 100.0 400.0 400.0 3
element elasticBeamColumn 47  7 21 20.0 29000.0 11154.0 100.0 400.0 400.0 1
element elasticBeamColumn 48 14 28 20.0 29000.0 11154.0 100.0 400.0 400.0 1

# Define gravity loads
# --------------------
# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {
    # Create nodal loads
    #    nd    FX  FY         FZ  MX  MY  MZ 
    load  2   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  9   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 16   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 23   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  3   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 10   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 17   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 24   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  4   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 11   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 18   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 25   0.0 0.0 [expr -$P] 0.0 0.0 0.0  
    load  5   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 12   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 19   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 26   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  6   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 13   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 20   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 27   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  7   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 14   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 21   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load 28   0.0 0.0 [expr -$P] 0.0 0.0 0.0
}
# ------------------------------
# End of model generation
# ------------------------------



# ------------------------------
# Start of analysis generation
# ------------------------------
# Create the system of equation
system BandGeneral
# Create the DOF numberer
numberer Plain
# Create the constraint handler
constraints Transformation
# Create the convergence test
test NormDispIncr 1.0e-12 10
# Create the integration scheme
integrator LoadControl 0.1
# Create the solution algorithm
algorithm Newton
# Create the analysis object
analysis Static
# ------------------------------
# End of analysis generation
# ------------------------------



# ------------------------------
# Start of recorder generation
# ------------------------------
# create a Recorder object for the nodal displacements at node 2
recorder Node -file Gravity_Dsp.out -time -nodeRange 1 28 -dof 1 2 3 4 5 6 disp
recorder Element -file Gravity_Frc.out -time -ele 1 2 3 4 force
# --------------------------------
# End of recorder generation
# --------------------------------



# ------------------------------
# Perform the gravity analysis
# ------------------------------
# perform the gravity load analysis, requires 10 steps to reach the load level
analyze 10
puts "\nGravity load analysis completed";

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0
remove recorders



# --------------------------------
# Perform an eigenvalue analysis
# --------------------------------
set pi [expr acos(-1.0)]
set lambda [eigen -fullGenLapack 42]
set omega1 [expr pow([lindex $lambda 0],0.5)]
puts "\nEigenvalues at start of transient:"
puts "|   lambda   |  omega   |  period | frequency |"
foreach lambda $lambda {
    set omega [expr pow($lambda,0.5)]
    set period [expr 2.0*$pi/$omega]
    set frequ [expr 1.0/$period]
    puts [format "| %5.3e | %8.4f | %7.4f | %9.4f |" $lambda $omega $period $frequ]
}



# ------------------------------
# Start of model generation
# ------------------------------

# Define dynamic loads
# --------------------
# Set time series to be passed to uniform excitation
set dt 0.005
set scale 0.75; # max = 0.85
set npts 9000
timeSeries Path 2 -filePath SCS052.txt -dt $dt -factor [expr $g*$scale]
timeSeries Path 3 -filePath SCS142.txt -dt $dt -factor [expr $g*$scale]
timeSeries Path 4 -filePath SCSUP.txt -dt $dt -factor [expr $g*$scale]

# Create UniformExcitation load pattern
#                         tag dir -accel tsTag 
pattern UniformExcitation  2   1  -accel   2
pattern UniformExcitation  3   2  -accel   3
pattern UniformExcitation  4   3  -accel   4

# calculate the Rayleigh damping factors for nodes & elements
set zeta 0.01
set beta [expr 2.0*$zeta/$omega1];
set alphaM     0.0;     # mass proportional damping;       D = alphaM*M
set betaK      0.0;     # stiffness proportional damping;  D = betaK*Kcurrent
set betaKinit  0.0;     # stiffness proportional damping;  D = beatKinit*Kinit
set betaKcomm  $beta;   # stiffness proportional damping;  D = betaKcomm*KlastCommit

# set the Rayleigh damping 
rayleigh $alphaM $betaK $betaKinit $betaKcomm
# ------------------------------
# End of model generation
# ------------------------------



# ------------------------------
# Start of recorder generation
# ------------------------------
# create a Recorder object for the nodal displacements at node 2
recorder Node -file Node_Dsp.out -time -node 2 3 4 5 6 7 29 30 31 32 33 34 -dof 1 2 3 4 5 6 disp
recorder Node -file Node_Vel.out -time -node 2 3 4 5 6 7 29 30 31 32 33 34 -dof 1 2 3 4 5 6 vel
recorder Node -file Node_Acc.out -time -node 2 3 4 5 6 7 29 30 31 32 33 34 -dof 1 2 3 4 5 6 accel
recorder Node -file Node_AbsAcc.out -timeSeries 2 3 4 -time -node 1 2 3 4 5 6 7 29 30 31 32 33 34 -dof 1 2 3 accel

recorder Element -file Elmt_Frc.out -time -ele 1 force
recorder Element -file Elmt_Def.out -time -ele 1 basicDeformation
recorder Element -file Elmt_N.out -time -ele 1 frictionModel normalForce
recorder Element -file Elmt_Vel.out -time -ele 1 frictionModel vel
recorder Element -file Elmt_Ff.out -time -ele 1 frictionModel frictionForce
recorder Element -file Elmt_COF.out -time -ele 1 frictionModel COF

# # recorder display "Display" xLoc yLoc xPixels yPixels -wipe 
recorder  display  "Display"  5  5  630  630 -wipe
# "normal" vector to the view window
# vpn-5.272000E-001  -6.871000E-001  +5.000000E-001
# "up" vector of the view window
# vup+0.000000E+000  +0.000000E+000  +1.000000E+000
# Projection Reference Point (direction vector to the eye)
# prp-4.965000E-001  -6.576000E-001  +5.666000E-001
# dimension of the view window
# viewWindow-6.390000E+002  +6.390000E+002  -6.390000E+002  +6.390000E+002
# center of the view window
# vrp+7.200000E+001  +7.200000E+001  +3.650000E+002
# display    elemDispOpt    nodeDispOpt    magFactor
# display1  3  +2.000000E+000
# --------------------------------
# End of recorder generation
# --------------------------------



# ------------------------------
# Start of analysis generation
# ------------------------------
# create the system of equation
system BandGeneral
# create the DOF numberer
numberer Plain
# create the constraint handler
constraints Transformation

# set the test parameters
set testType NormDispIncr
set testTol 1.0e-12;
set testIter 50;
test $testType $testTol $testIter

# set the integrator parameters
integrator Newmark 0.5 0.25

# set the algorithm parameters
set algoType KrylovNewton
algorithm  $algoType

# create the analysis object 
analysis Transient
# ------------------------------
# End of analysis generation
# ------------------------------



# ------------------------------
# Finally perform the analysis
# ------------------------------
logFile "TestFPS3d_4.log"

set dtAna [expr $dt/2.5]
set dtMin 1.0e-8
set dtMax $dtAna

set ok 0;
set tFinal [expr $npts * $dt]
set tCurrent [getTime "%1.12E"]

record
while {$ok == 0 && $tCurrent < $tFinal} {
    
    set ok [analyze 1 $dtAna]
    
    if {$ok != 0} {
        if {[expr $dtAna/2.0] >= $dtMin} {
            set dtAna [expr $dtAna/2.0]
            puts [format "\nREDUCING time step size (dtNew = %1.6e)" $dtAna]
            set ok 0
        }
    } else {
        set tCurrent [getTime "%1.12E"]
        puts [format "t = %1.4f sec" $tCurrent]
        if {[expr $dtAna*2.0] <= $dtMax} {
            set dtAna [expr $dtAna*2.0]
            puts [format "\nINCREASING time step size (dtNew = %1.6e)" $dtAna]
        }
    }
}

if {$ok != 0} {
    puts [format "\nModel failed (time = %1.3e)" $tCurrent]
} else {
    puts [format "\nResponse-history analysis completed"]
}

wipe
exit
# --------------------------------
# End of analysis
# --------------------------------
