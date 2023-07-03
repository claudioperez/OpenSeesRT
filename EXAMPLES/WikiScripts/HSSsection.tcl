############################################################################################
# HSSsection.tcl
#
# This routine creates a fiber section: AISC standard HSS section 
# 
# Variables
	# secID = section ID number
	# matID = material ID number 
	# d  = nominal depth	
	# t  = tube tickness
	# nfdy = number of fibers along depth that goes along local y axis 
	# nfty = number of fibers along thickness that goes along local y axis
	# nfdz = number of fibers along depth that goes along local z axis
	# nftz = number of fibers along thickness that goes along local z axis
############################################################################################

############################################################################################

proc HSSsection { secID matID d t nfdy nfty nfdz nftz} {
	set dw [expr $d - 2 * $t]
	set y1 [expr -$d/2]
	set y2 [expr -$dw/2]
	set y3 [expr  $dw/2]
	set y4 [expr  $d/2]
  
	set z1 [expr -$d/2]
	set z2 [expr -$dw/2]
	set z3 [expr  $dw/2]
	set z4 [expr  $d/2]
  
	section fiberSec  $secID  {
   		#                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
   		patch quadr  $matID  $nftz $nfdy   $y2 $z4   $y2 $z3   $y3 $z3   $y3 $z4
   		patch quadr  $matID  $nftz $nfdy   $y2 $z2   $y2 $z1   $y3 $z1   $y3 $z2
   		patch quadr  $matID  $nfdz $nfty   $y1 $z4   $y1 $z1   $y2 $z1   $y2 $z4
   		patch quadr  $matID  $nfdz $nfty   $y3 $z4   $y3 $z1   $y4 $z1   $y4 $z4
	}
}