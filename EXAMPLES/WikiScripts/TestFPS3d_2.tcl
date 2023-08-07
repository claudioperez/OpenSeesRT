# File: TestFPS3d_2.tcl
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
# singleFPBearing element. It models an isolated one story
# stick and the bearing element has finite length. It also
# tests the different friction models.

# ------------------------------
# Start of model generation
# ------------------------------
# Create ModelBuilder (with three-dimensions and 6 DOF/node)
model BasicBuilder -ndm 3 -ndf 6

# Define geometry for model
# -------------------------
set g [expr 32.174*12.0]
set P 9.0
set mass [expr $P/$g]
#    tag   xCrd   yCrd   zCrd        mass
node  1     0.0    0.0    0.0
node  2     0.0    0.0   10.0  -mass $mass $mass $mass 0.0 0.0 0.0
node  3     0.0    0.0  154.0  -mass $mass $mass $mass 0.0 0.0 0.0

# Set the boundary conditions
#   tag DX DY DZ RX RY RZ
fix 1   1  1  1  1  1  1
fix 2   0  0  0  1  1  1
fix 3   0  0  0  1  1  1

# Define material models
# ----------------------
set mv [expr 2.0*$mass]
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
geomTransf Linear  1  0 1 0

# Define elements
# ---------------
# element flatSliderBearing eleTag NodeI NodeJ frnMdlTag kInit -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element flatSliderBearing 1 1 2 1 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0

# element singleFPBearing eleTag NodeI NodeJ frnMdlTag Reff kInit -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
element singleFPBearing 1 1 2 1 34.68 250.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0

# element RJWatsonEqsBearing eleTag NodeI NodeJ frnMdlTag kInit k2 k3 mu -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element RJWatsonEqsBearing 1 1 2 1 250.0 0.519 0.0 3.0 -P 1 -T 2 -My 2 -Mz 2 -orient 1 0 0

# element elasticBeamColumn eleTag NodeI NodeJ A E G J Iy Iz geoTranTag <-mass massDens> 
element elasticBeamColumn 2 2 3 20.0 29000.0 11154.0 100.0 200.0 200.0 1 

# Define gravity loads
# --------------------
# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {
    # Create nodal loads
    #    nd    FX  FY         FZ  MX  MY  MZ 
    load  2   0.0 0.0 [expr -$P] 0.0 0.0 0.0
    load  3   0.0 0.0 [expr -$P] 0.0 0.0 0.0
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
constraints Plain
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
recorder Node -file Gravity_Dsp.out -time -node 2 3 -dof 1 2 3 4 5 6 disp
recorder Element -file Gravity_Frc.out -time -ele 1 force
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
set lambda [eigen -fullGenLapack 3]
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
set scale 1.0; # max = 1.1
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
recorder Node -file Node_Dsp.out -time -node 2 3 -dof 1 2 3 4 5 6 disp
recorder Node -file Node_Vel.out -time -node 2 3 -dof 1 2 3 4 5 6 vel
recorder Node -file Node_Acc.out -time -node 2 3 -dof 1 2 3 4 5 6 accel
recorder Node -file Node_AbsAcc.out -timeSeries 2 3 4 -time -node 1 2 3 -dof 1 2 3 accel

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
# viewWindow-1.350000E+002  +1.350000E+002  -1.350000E+002  +1.350000E+002
# center of the view window
# vrp+0.000000E+000  +0.000000E+000  +7.700000E+001
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
constraints Plain

# set the test parameters
set testType NormDispIncr
set testTol 1.0e-12;
set testIter 25;
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
logFile "TestFPS3d_2.log"

set dtAna [expr $dt/1.0]
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
