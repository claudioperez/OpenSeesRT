###########################################################################################################
# rotLeaningCol.tcl
# Procedure that creates a zero-stiffness elastic rotational spring for the leaning column
# and constrains the translations DOFs of the spring.
#
# Written by: Laura Eads
# Date: 07/16/2010
#
# Formal arguments
#       eleID   - unique element ID for this zero length rotational spring
#       nodeR   - node ID which will be retained by the multi-point constraint
#       nodeC   - node ID which will be constrained by the multi-point constraint
#
##########################################################################################################

proc rotLeaningCol {eleID nodeR nodeC} {

	#Spring Stiffness
	set K 1e-9; # k/in

	# Create the material and zero length element (spring)
    uniaxialMaterial Elastic  $eleID  $K	
	element zeroLength $eleID $nodeR $nodeC -mat $eleID -dir 6

	# Constrain the translational DOF with a multi-point constraint	
	#   		retained constrained DOF_1 DOF_2 
	equalDOF    $nodeR     $nodeC     1     2
}