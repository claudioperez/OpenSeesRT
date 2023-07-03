############################################################################################
# rotSpring2DModIKModel.tcl
#
# This routine creates a uniaxial material spring with deterioration
# 
# Spring follows: Bilinear Response based on Modified Ibarra Krawinkler Deterioration Model 
#
# Written by: Dimitrios G. Lignos, Ph.D.
#
# Variables
# 	$eleID = 	Element Identification (integer)
# 	$nodeR = 	Retained/master node
# 	$nodeC = 	Constrained/slave node
# 	$K = 		Initial stiffness after the modification for n (see Ibarra and Krawinkler, 2005)
# 	$asPos = 	Strain hardening ratio after n modification (see Ibarra and Krawinkler, 2005)
# 	$asNeg = 	Strain hardening ratio after n modification (see Ibarra and Krawinkler, 2005)
# 	$MyPos = 	Positive yield moment (with sign)
# 	$MyNeg = 	Negative yield moment (with sign)
# 	$LS = 		Basic strength deterioration parameter (see Lignos and Krawinkler, 2009)
# 	$LK = 		Unloading stiffness deterioration parameter (see Lignos and Krawinkler, 2009)
# 	$LA = 		Accelerated reloading stiffness deterioration parameter (see Lignos and Krawinkler, 2009)
# 	$LD = 		Post-capping strength deterioration parameter (see Lignos and Krawinkler, 2009)
# 	$cS = 		Exponent for basic strength deterioration
# 	$cK = 		Exponent for unloading stiffness deterioration
# 	$cA = 		Exponent for accelerated reloading stiffness deterioration
# 	$cD = 		Exponent for post-capping strength deterioration
# 	$th_pP = 	Plastic rotation capacity for positive loading direction
# 	$th_pN = 	Plastic rotation capacity for negative loading direction
# 	$th_pcP = 	Post-capping rotation capacity for positive loading direction
# 	$th_pcN = 	Post-capping rotation capacity for negative loading direction
# 	$ResP = 	Residual strength ratio for positive loading direction
# 	$ResN = 	Residual strength ratio for negative loading direction
# 	$th_uP = 	Ultimate rotation capacity for positive loading direction
# 	$th_uN = 	Ultimate rotation capacity for negative loading direction
# 	$DP = 		Rate of cyclic deterioration for positive loading direction
# 	$DN = 		Rate of cyclic deterioration for negative loading direction
#
# References:
#		Ibarra, L. F., and Krawinkler, H. (2005). “Global collapse of frame structures under seismic excitations,” Technical Report 152, The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
# 		Ibarra, L. F., Medina, R. A., and Krawinkler, H. (2005). “Hysteretic models that incorporate strength and stiffness deterioration,” International Journal for Earthquake Engineering and Structural Dynamics, Vol. 34, No.12, pp. 1489-1511.
# 		Lignos, D. G., and Krawinkler, H. (2010). “Deterioration Modeling of Steel Beams and Columns in Support to Collapse Prediction of Steel Moment Frames”, ASCE, Journal of Structural Engineering (under review).
# 		Lignos, D. G., and Krawinkler, H. (2009). “Sidesway Collapse of Deteriorating Structural Systems under Seismic Excitations,” Technical Report 172, The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
#
############################################################################################
#
proc rotSpring2DModIKModel {eleID nodeR nodeC K asPos asNeg MyPos MyNeg LS LK LA LD cS cK cA cD th_pP th_pN th_pcP th_pcN ResP ResN th_uP th_uN DP DN} {
#
# Create the zero length element
      uniaxialMaterial Bilin  $eleID  $K  $asPos $asNeg $MyPos $MyNeg $LS $LK $LA $LD $cS $cK $cA $cD $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN $DP $DN;

      element zeroLength $eleID $nodeR $nodeC -mat $eleID -dir 6

# Constrain the translational DOF with a multi-point constraint
#          			retained constrained DOF_1 DOF_2 ... DOF_n
           equalDOF    $nodeR     $nodeC     1     2
}