# --------------------------------------------------------------------------------------------------
# Example 5. 2D Frame --  Static Reversed Cyclic Analysis
#		Silvia Mazzoni & Frank McKenna, 2006
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
pattern Plain 200 Linear {;			# define load pattern
	for {set level 2} {$level <=[expr $NStory+1]} {incr level 1} {
		set Fi [lindex $iFi [expr $level-1-1]];		# lateral load coefficient
		for {set pier 1} {$pier <= [expr $NBay+1]} {incr pier 1} {
			set nodeID [expr $level*10+$pier]
			load $nodeID $Fi 0.0 0.0 0.0 0.0 0.0
		}
	}
};		# end load pattern

# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

#  ---------------------------------    perform Static Cyclic Displacements Analysis
set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis
set count 0
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
			incr count
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
