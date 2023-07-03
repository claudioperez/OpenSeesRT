# Example analysis file for LimitStateMaterial
# 
# Units: kip, in
# KJE, Feb 2003
# modified by Mohammad Reza Azadi Kakavand to accomodate it for transient analysis, Nov 2011

wipe

source tags.tcl
source DisplayModel2D.tcl
source DisplayPlane.tcl

###########################
# Build model
###########################

model BasicBuilder -ndm 2 -ndf 3
puts "HI"
set dataDir Analysis_LimitStateMaterial(Disp)-Dynamic

file mkdir $dataDir
################################
# Define nodal mesh and B.C.s
################################
set L 58.0

#    tag  X   Y
node  1  0.0 0.0
node  2  0.0 0.0
node  3  0.0 $L
node  4  0.0  $L
node  5  0.0  $L


#   tag DX DY RZ
fix  1   1  1  1
fix  4   0  0  1
fix  5   1  1  1

mass 4 0.1813 1.0e-9 0
##############################
# Create column section
##############################

source CenterColSecFiber.tcl


#################################
# Define the beam-column element
#################################
geomTransf PDelta 1
set nint 5
element nonlinearBeamColumn $bcTag 2 3 $nint $flexSec 1 -iter 5 1e-15


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
recorder Node -file $dataDir/RX.out -time -node 1 -dof 1 reaction

recorder Node -file $dataDir/RY.out -time -node 1 -dof 2 reaction

recorder Drift -file $dataDir/Drift.out -time -iNode 2  -jNode 3 -dof 1 -perpDirn 2;		# lateral drift


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
#system ProfileSPD

system BandGeneral
#constraints Penalty 1.0e12 1.0e12

constraints Plain
numberer Plain
test NormDispIncr 1.0e-3 25 0
algorithm Newton
integrator LoadControl 0 1 0 0
analysis Static

# Apply gravity load in one step
analyze 1
wipeAnalysis

source Dynamic.EQ.Uniform_LimitState.tcl


