############################################################################################
# RectangularRCsection2D.tcl
#
# This routine creates a fiber section: 
# 
# Variables
	# secID = section ID number
	# matID = material ID number 
	# d  = section depth	
	# b = section width	
	# cover = concrete cover
	# nfCoreY = number of fibers for concrete in y-direction -- core concrete
	# nfCoreZ = number of fibers for concrete in z-direction -- core concrete
	# nfCoverYlong  = number of fibers for long cover patch of cover in y-direction
	# nfCoverYshort = number of fibers for short cover patch of cover in y-direction
	# As = area of one longitudinal bar
	# nBarsLayer1 = number of bars in the top layer
	# nBarsLayer2 = number of bars in the middle layer
	# nBarsLayer3 = number of bars in the bottom layer
############################################################################################

############################################################################################

proc RectangularRCsection2D { secID IDcore IDcover IDreinf d b cover nfCoreY nfCoreZ nfCoverYlong nfCoverYshort As nBarsLayer1 nBarsLayer2 nBarsLayer3} {

	set y1 [expr $d/2.0]
	set z1 [expr $b/2.0]

	section Fiber $secID {

		# Create the concrete core fibers
		# 			
		patch rect $IDcore $nfCoreY $nfCoreZ [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

		# Create the concrete cover fibers (left, right, bottom, top)                                                                                                  
		patch rect $IDcover $nfCoverYlong  1  [expr -$y1] [expr $z1-$cover] $y1 $z1
		patch rect $IDcover $nfCoverYlong  1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
		patch rect $IDcover $nfCoverYshort 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
		patch rect $IDcover $nfCoverYshort 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

		# Create the reinforcing fibers (top, middle, bottom)                                                                                                          
		layer straight $IDreinf $nBarsLayer1 $As [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
		layer straight $IDreinf $nBarsLayer2 $As 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
		layer straight $IDreinf $nBarsLayer3 $As [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]
	}	
  
}