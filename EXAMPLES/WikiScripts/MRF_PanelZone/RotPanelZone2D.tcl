###########################################################################################################
# rotPanelZone2D.tcl
# Procedure that creates a rotational spring and constrains the corner nodes of a panel zone
# 
# The equations and process are based on: Krawinkler Model for Panel Zones
# Reference:  Gupta, A., and Krawinkler, H. (1999). "Seismic Demands for Performance Evaluation of Steel Moment Resisting Frame Structures,"
#            Technical Report 132, The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
#
#
# Written by: Dimitrios Lignos
# Date: 11/09/2008
#
# Formal arguments
#       eleID   - unique element ID for this zero length rotational spring
#       nodeR   - node ID which will be retained by the multi-point constraint
#       nodeC   - node ID which will be constrained by the multi-point constraint
#       E       - modulus of elasticity
#       Fy      - yield strength
#       dc      - column depth
#       bf_c    - column flange width
#       tf_c    - column flange thickness
#       tp      - panel zone thickness
#       db      - beam depth
#       Ry      - expected value for yield strength --> Typical value is 1.2
#       as      - assumed strain hardening
##########################################################################################################

proc rotPanelZone2D {eleID nodeR nodeC E Fy dc bf_c tf_c tp db Ry as} {

# Trilinear Spring
# Yield Shear
	set Vy [expr 0.55 * $Fy * $dc * $tp];
# Shear Modulus
	set G [expr $E/(2.0 * (1.0 + 0.30))]
# Elastic Stiffness
	set Ke [expr 0.95 * $G * $tp * $dc];
# Plastic Stiffness
	set Kp [expr 0.95 * $G * $bf_c * ($tf_c * $tf_c) / $db];

# Define Trilinear Equivalent Rotational Spring
# Yield point for Trilinear Spring at gamma1_y
	set gamma1_y [expr $Vy/$Ke]; set M1y [expr $gamma1_y * ($Ke * $db)];
# Second Point for Trilinear Spring at 4 * gamma1_y
	set gamma2_y [expr 4.0 * $gamma1_y]; set M2y [expr $M1y + ($Kp * $db) * ($gamma2_y - $gamma1_y)];
# Third Point for Trilinear Spring at 100 * gamma1_y
	set gamma3_y [expr 100.0 * $gamma1_y]; set M3y [expr $M2y + ($as * $Ke * $db) * ($gamma3_y - $gamma2_y)];
  
  
# Hysteretic Material without pinching and damage (same mat ID as Ele ID)
    uniaxialMaterial Hysteretic $eleID $M1y $gamma1_y  $M2y $gamma2_y $M3y $gamma3_y [expr -$M1y] [expr -$gamma1_y] [expr -$M2y] [expr -$gamma2_y] [expr -$M3y] [expr -$gamma3_y] 1 1 0.0 0.0 0.0
	
	element zeroLength $eleID $nodeR $nodeC -mat $eleID -dir 6

	equalDOF    $nodeR     $nodeC     1     2
	# Constrain the translational DOF with a multi-point constraint
	# Left Top Corner of PZ
	set nodeR_1 [expr $nodeR - 2];
	set nodeR_2 [expr $nodeR_1 + 1];
	# Right Bottom Corner of PZ
	set nodeR_6 [expr $nodeR + 3];
	set nodeR_7 [expr $nodeR_6 + 1];
	# Left Bottom Corner of PZ
	set nodeL_8 [expr $nodeR + 5];
	set nodeL_9 [expr $nodeL_8 + 1];
	#          retained constrained DOF_1 DOF_2 
	equalDOF    $nodeR_1     $nodeR_2    1     2
	equalDOF    $nodeR_6     $nodeR_7    1     2
	equalDOF    $nodeL_8     $nodeL_9    1     2
}