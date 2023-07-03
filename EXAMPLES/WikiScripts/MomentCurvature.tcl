proc MomentCurvature {secTag axialLoad maxK {numIncr 100} } {
	# Define two nodes at (0,0)
	node 1 0.0 0.0
	node 2 0.0 0.0

	# Fix all degrees of freedom except axial and bending
	fix 1 1 1 1
	fix 2 0 1 0

	# Define element
	#                         tag ndI ndJ  secTag
	element zeroLengthSection  1   1   2  $secTag

	# Create recorder
	recorder Node -file section$secTag.out -time -node 2 -dof 3 disp

	# Define constant axial load
	pattern Plain 1 "Constant" {
		load 2 $axialLoad 0.0 0.0
	}

	# Define analysis parameters
	integrator LoadControl 0.0
	system SparseGeneral -piv;	# Overkill, but may need the pivoting!
	test NormUnbalance 1.0e-9 10
	numberer Plain
	constraints Plain
	algorithm Newton
	analysis Static

	# Do one analysis for constant axial load
	analyze 1

	# Define reference moment
	pattern Plain 2 "Linear" {
		load 2 0.0 0.0 1.0
	}

	# Compute curvature increment
	set dK [expr $maxK/$numIncr]

	# Use displacement control at node 2 for section analysis
	integrator DisplacementControl 2 3 $dK 1 $dK $dK

	# Do the section analysis
	analyze $numIncr
}


# Define model builder
# --------------------
model basic -ndm 2 -ndf 3

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014

# Cover concrete (unconfined)
uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006

# STEEL
# Reinforcing steel 
set fy 60.0;      # Yield stress
set E 30000.0;    # Young's modulus
#                        tag  fy E0    b
uniaxialMaterial Steel01  3  $fy $E 0.01

# Define cross-section for nonlinear columns
# ------------------------------------------

# set some paramaters
set colWidth 15
set colDepth 24 

set cover  1.5
set As    0.60;     # area of no. 7 bars

# some variables derived from the parameters
set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber 1 {

    # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect 2 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect 2  2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect 2  2 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $As [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
    layer straight 3 2 $As 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
    layer straight 3 3 $As [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]

}    

# Estimate yield curvature
# (Assuming no axial load and only top and bottom steel)
set d [expr $colDepth-$cover]	;# d -- from cover to rebar
set epsy [expr $fy/$E]	;# steel yield strain
set Ky [expr $epsy/(0.7*$d)]

# Print estimate to standard output
puts "Estimated yield curvature: $Ky"

# Set axial load 
set P -180

set mu 15;		# Target ductility for analysis
set numIncr 100;	# Number of analysis increments

# Call the section analysis procedure
MomentCurvature 1 $P [expr $Ky*$mu] $numIncr