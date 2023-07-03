# Example analysis file for LimitStateMaterial
# 
# Units: kip, in
# KJE, Feb 2003
# debugged by Mohammad Reza Azadi Kakavand, Nov 2011

source Tags.tcl

# Number of analysis steps and step size
#set nSteps 35000
set nSteps 6000
set dlamda 0.1

###########################
# Build model
###########################

model BasicBuilder -ndm 2 -ndf 3
puts "HI"
set dataDir Analysis_LimitStateMaterial(Load)-Pushover
file mkdir $dataDir
################################
# Define nodal mesh and B.C.s
################################
set L 58.0

#    tag  X   Y
node  1  0.0 0.0
node  2  0.0 0.0
node  3  0.0  $L
node  4  0.0  $L
node  5  0.0  $L

#   tag DX DY RZ
fix  1   1  1  1
fix  4   0  0  1
fix  5   1  1  1


##############################
# Create column section
##############################

source CenterColSecFiber.tcl


#################################
# Define the beam-column element
#################################
geomTransf PDelta 1
set nint 5
element forceBeamColumn $bcTag 2 3 $nint $flexSec 1 -iter 5 1e-15


####################################
# Define the zero-length end springs
####################################
#rigid material
uniaxialMaterial Elastic $rigidMatTag 9.9e9

#bottom of column slip spring 
element zeroLength 1 1 2 -mat $rigidMatTag $rigidMatTag $centerSlipTag -dir 1 2 6
	
#top of column springs incl shear and axial
element zeroLength 3 3 4 -mat $shearTag $axialFailTag $centerSlipTag -dir 1 2 6 

#soft spring to take axial load after failure
uniaxialMaterial Elastic $softMatTag 99.9
element zeroLength 4 4 5 -mat $softMatTag -dir 2


#####################
# Recorders
#####################

# record displacements for node with applied displacement
recorder Node -file $dataDir/node4DispX.out -time -node 4 -dof 1 disp
recorder Node -file $dataDir/node4DispY.out -time -node 4 -dof 2 disp
recorder Node -file $dataDir/RX.out -time -node 5 -dof 1 reaction
recorder Node -file $dataDir/RY.out -time -node 5 -dof 2 reaction
# record end section forces and deformations
recorder Element -file $dataDir/secforce1Column.out -time  -ele $bcTag  section 1 force
recorder Element -file $dataDir/secdeform1Column.out -time  -ele $bcTag section 1 deformation

# record beam-column element forces and rotations
recorder Element -file $dataDir/axialSpr.out  -time  -ele 4 force
recorder Element -file $dataDir/eleforcebasic.out -time -ele $bcTag force
recorder Element -file $dataDir/elerotbasic.out -time  -ele $bcTag rotation


###########################
# Apply gravity loads
###########################

# Initial axial load 
set P -70.0 

# Constant gravity load
pattern Plain 1 Constant {
	load 4 0.0 $P 0.0 
	}

# Analysis options for gravity loads
initialize
system ProfileSPD
constraints Penalty 1.0e12 1.0e12
numberer RCM
test NormDispIncr 1.0e-8 25 0
algorithm Newton
integrator LoadControl 0 1 0 0
analysis Static

# Apply gravity load in one step
analyze 1
wipeAnalysis

###### Loading option ######
set loading push
#set loading cyclic
#
# Analysis control option
set control load
#set control displ
#################################
# Apply transverse displacements
#################################
if {$loading == "push"} {

    # Define displacements of node 4 for pushover
    set du 0.01
    pattern Plain 2 "Linear -factor $du" {
	sp 4 1 1.0 
    }
	
} elseif {$loading == "cyclic"} {

    set path {1.2557524028480316e-003  1.2557524028480316e-003  3.0608252487525078e-003  1.5641001656305775e-003  1.1015791633965932e-003 -1.8313432524451514e+000}

    pattern Plain 2 "Series -dt 1.0 -values $path" {
		sp 4 1 1.0 }
	
} else {
	puts stderr "Invalid loading option: $loading"
}


###############################################
# Set Analysis options for transverse loading
###############################################

# Create the system of equation
system BandGeneral

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the integration scheme, 
if {$control == "load"} {
	#LoadControl scheme using constant steps of dlamda
	#                      dlamda1 Jd minLamda maxLamda
	integrator LoadControl $dlamda 1  $dlamda $dlamda
	constraints Penalty 1.0e14  1.0e14
	
} elseif {$control == "displ"} {
        #DisplacementControl scheme with constant loading rate
        #                              node dof du1 Jd mindu maxdu
      
      integrator DisplacementControl 2     1  $du 1  $du   $du
	constraints Penalty 1.0e14  1.0e14

} else {
	puts stderr "Invalid control option: $control"
}
	
# Create the convergence test, Norm of the displacement increment
#                 tol     maxNumIter printFlag
test NormDispIncr 1.0e-10 25         0

# Create the solution algorithm, a Newton-Rahpson algorithm is created
algorithm Newton
#algorithm KrylovNewton

# create the analysis object 
analysis Static 


#########################################
# Analyze model with transverse loading
#########################################
set ok 0
set n 1
while {$n < $nSteps && $ok == 0} {
	set ok [analyze 1]
	if {$ok != 0} {
		test NormDispIncr 1.0e-10 10000 1     ;# increase max number of iterations 
		algorithm ModifiedNewton -initial      ;# use initial elastic stiffness for NewtonRaphson
		#puts "Time Newton Initial [getTime]" ;# output time algorithm was changed
		set ok [analyze 1 ] 			;# analyse 1 step with new settings
		#set ans [gets stdin] 		        ;# pauses tcl script
		algorithm Newton                        ;# restore algorithm to Newton with current stiffness
		test NormDispIncr 1.0e-10 25 0          ;# restore max number of iterations
	}
	incr n
}

if {$ok != 0} {
	puts "ANALYSIS FAILED"
} else {
	puts "SUCCESSFULL"
}

wipe
