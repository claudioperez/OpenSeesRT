# File: TestSlider2d_3.tcl
#
# $Revision: $
# $Date: $
# $URL: $
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Created: 02/09
# Revision: A
#
# Purpose: this file tests the 2D flatSliderBearing or the
# singleFPBearing element. It models an isolated one story  
# one bay building and the bearing element has finite length.
# It also tests the different friction models.

# ------------------------------
# Start of model generation
# ------------------------------
# Create ModelBuilder (with two-dimensions and 3 DOF/node)
model BasicBuilder -ndm 2 -ndf 3

# Define geometry for model
# -------------------------
set g [expr 32.174*12.0]
set P 9.0
set mass [expr $P/$g]
#    tag   xCrd   yCrd        mass
node  1     0.0    0.0
node  2     0.0   10.0  -mass $mass $mass 0.0
node  3     0.0  154.0  -mass $mass $mass 0.0
node  4   144.0    0.0
node  5   144.0   10.0  -mass $mass $mass 0.0
node  6   144.0  154.0  -mass $mass $mass 0.0

# set the boundary conditions
#   tag DX DY RZ
fix  1   1  1  1
fix  4   1  1  1

# Define material models
# ----------------------
set mv [expr 2.0*$mass] 
set kv 7500.0
set zetaVertical 0.02
set cv [expr 2.0*$zetaVertical*sqrt($kv*$mv)]
uniaxialMaterial Elastic 1 $kv $cv
#uniaxialMaterial ElasticMultiLinear 1 $cv -strain -0.002 -0.001 -0.000875 -0.00075 -0.000625 -0.0005 -0.000375 -0.00025 -0.000125 0.000125 -stress -12.5 -5 -4.0924 -3.2476 -2.4705 -1.7678 -1.1482 -0.625 -0.22097 0.22097
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
# geomTransf type tag 
geomTransf Linear 1 

# Define elements
# ---------------
# element flatSliderBearing eleTag NodeI NodeJ frnMdlTag kInit -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
element flatSliderBearing 1 1 2 1 250.0 -P 1 -Mz 2
element flatSliderBearing 2 4 5 1 250.0 -P 1 -Mz 2

# element singleFPBearing eleTag NodeI NodeJ frnMdlTag Reff kInit -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element singleFPBearing 1 1 2 1 34.68 250.0 -P 1 -Mz 2
#element singleFPBearing 2 4 5 1 34.68 250.0 -P 1 -Mz 2

# element RJWatsonEqsBearing eleTag NodeI NodeJ frnMdlTag kInit k2 k3 mu -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>
#element RJWatsonEqsBearing 1 1 2 1 250.0 0.519 0.0 3.0 -P 1 -Mz 2
#element RJWatsonEqsBearing 2 4 5 1 250.0 0.519 0.0 3.0 -P 1 -Mz 2

# element elasticBeamColumn eleTag NodeI NodeJ A E Iz geoTranTag <alpha d> <-mass massDens> 
element elasticBeamColumn 3 2 3 20.0 29000.0 400.0 1 
element elasticBeamColumn 4 5 6 20.0 29000.0 400.0 1 
element elasticBeamColumn 5 2 5 20.0 29000.0 400.0 1 
element elasticBeamColumn 6 3 6 20.0 29000.0 400.0 1 

# Define gravity loads
# --------------------
# create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {
    # create nodal loads
    #    nd    FX          FY   MZ 
    load  2   0.0  [expr -$P]  0.0
    load  3   0.0  [expr -$P]  0.0
    load  5   0.0  [expr -$P]  0.0
    load  6   0.0  [expr -$P]  0.0
}
# ------------------------------
# End of model generation
# ------------------------------



# ------------------------------
# Start of analysis generation
# ------------------------------
# create the system of equation
system BandGeneral
# create the DOF numberer
numberer Plain
# create the constraint handler
constraints Plain
# create the convergence test
test NormDispIncr 1.0e-12 10
# create the integration scheme
integrator LoadControl 0.1
# create the solution algorithm
algorithm Newton
# create the analysis object
analysis Static
# ------------------------------
# End of analysis generation
# ------------------------------



# ------------------------------
# Start of recorder generation
# ------------------------------
# create a Recorder object for the nodal displacements at node 2
recorder Node -file Gravity_Dsp.out -time -node 2 3 5 6 -dof 1 2 3 disp
recorder Element -file Gravity_Frc.out -time -ele 1 2 force
# --------------------------------
# End of recorder generation
# --------------------------------



# ------------------------------
# Perform the gravity analysis
# ------------------------------
# perform the gravity load analysis, requires 10 steps to reach the load level
analyze 10
puts "\nGravity load analysis completed";

# set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0
remove recorders



# --------------------------------
# Perform an eigenvalue analysis
# --------------------------------
set pi [expr acos(-1.0)]
set lambda [eigen -fullGenLapack 8]
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
# set time series to be passed to uniform excitation
set dt 0.005
set scale 1.0; # max = 1.7
set npts 9000
timeSeries Path 2 -filePath SCS052.AT2 -dt $dt -factor [expr $g*$scale]
timeSeries Path 3 -filePath SCSUP.AT2 -dt $dt -factor [expr $g*$scale]

# create UniformExcitation load pattern
#                         tag dir -accel tsTag 
pattern UniformExcitation  2   1  -accel   2
pattern UniformExcitation  3   2  -accel   3

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
recorder Node -file Node_Dsp.out -time -node 2 3 5 6 -dof 1 2 3 disp
recorder Node -file Node_Vel.out -time -node 2 3 5 6 -dof 1 2 3 vel
recorder Node -file Node_Acc.out -time -node 2 3 5 6 -dof 1 2 3 accel
recorder Node -file Node_AbsAcc.out -timeSeries 2 3 -time -node 1 2 3 5 6 -dof 1 2 accel

recorder Element -file Elmt_Frc.out -time -ele 1 2 force
recorder Element -file Elmt_Def.out -time -ele 1 2 basicDeformation
recorder Element -file Elmt_N.out -time -ele 1 2 frictionModel normalForce
recorder Element -file Elmt_Vel.out -time -ele 1 2 frictionModel vel
recorder Element -file Elmt_Ff.out -time -ele 1 2 frictionModel frictionForce
recorder Element -file Elmt_COF.out -time -ele 1 2 frictionModel COF

# recorder display "Display" xLoc yLoc xPixels yPixels -wipe 
recorder  display  "Display"  5  5  630  630 -wipe
# "normal" vector to the view window
vpn +0.000000E+000  +0.000000E+000  +1.000000E+000
# "up" vector of the view window
vup +0.000000E+000  +1.000000E+000  +0.000000E+000
# Projection Reference Point (direction vector to the eye)
prp +0.000000E+000  +0.000000E+000  +1.000000E+000
# dimension of the view window
viewWindow -1.160000E+002  +1.160000E+002  -1.160000E+002  +1.160000E+002
# center of the view window
vrp +7.200000E+001  +7.700000E+001  +0.000000E+000
# display    elemDispOpt    nodeDispOpt    magFactor
display 1  3  +2.000000E+000
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
set algoType Newton
algorithm  $algoType

# create the analysis object 
analysis Transient
# ------------------------------
# End of analysis generation
# ------------------------------



# ------------------------------
# Finally perform the analysis
# ------------------------------
logFile "TestSlider2d_3.log"

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
