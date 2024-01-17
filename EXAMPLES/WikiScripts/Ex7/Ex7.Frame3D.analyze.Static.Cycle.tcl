# --------------------------------------------------------------------------------------------------
# # Example 7. Static Reversed Cycles
#                       Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# source in procedures
source LibGeneratePeaks.tcl;		# procedure to generate displacement increments for cyclic peaks

# characteristics of cyclic analysis	
set iDmax "0.005  0.01 0.025 0.05 0.1";	# vector of displacement-cycle peaks, in terms of storey drift ratio
set Fact $LBuilding ;			# scale drift ratio by storey height for displacement cycles
set Dincr [expr 0.001*$LBuilding ];	# displacement increment for pushover. you want this to be very small, but not too small to slow analysis
set CycleType Full;			# you can do Full / Push / Half cycles with the proc
set Ncycles 1;			# specify the number of cycles at each peak

# -- STATIC PUSHOVER/CYCLIC ANALYSIS
# create load pattern for lateral pushover load coefficient when using linear load pattern
# need to apply lateral load only to the master nodes of the rigid diaphragm at each floor
pattern Plain 200 Linear {;			
	load 1121 $F2 0. 0. 0. 0. 0.
	load 1131 $F3 0. 0. 0. 0. 0.
	load 1141 $F4 0. 0. 0. 0. 0.
}

# Define DISPLAY -------------------------------------------------------------
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 2;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels
#recorder plot $dataDir/DFree.out Displ-X [expr $xPixels+10] 10 300 300 -columns 2 1; # a window to plot the nodal displacements versus time

# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

#  ---------------------------------    perform Static Cyclic Displacements Analysis
set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis
foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType $Fact];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep
			set Dincr [expr $D1 - $D0]
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			analysis Static
			# ----------------------------------------------first analyze command------------------------
			set ok [analyze 1]
			# ----------------------------------------------if convergence failure-------------------------
			if {$ok != 0} {
				# if analysis fails, we try some other stuff
				# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
				if {$ok != 0} {
					puts "Trying Newton with Initial Tangent .."
					test NormDispIncr   $Tol 2000 0
					algorithm Newton -initial
					set ok [analyze 1]
					test $testTypeStatic $TolStatic      $maxNumIterStatic    0
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					puts "Trying Broyden .."
					algorithm Broyden 8
					set ok [analyze 1 ]
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					puts "Trying NewtonWithLineSearch .."
					algorithm NewtonLineSearch 0.8 
					set ok [analyze 1]
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					puts $putout
					return -1
				}; # end if
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}