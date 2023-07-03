# --------------------------------------------------------------------------------------------------
# Example4. 2D Portal Frame--  Static Reversed-Cyclic Analysis
#                                Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# we need to set up parameters that are particular to the model.
set IDctrlNode 3;				# node where displacement is read for displacement control
set IDctrlDOF 1;				# degree of freedom of displacement read for displacement contro
# characteristics of cyclic analysis	
set iDmax "0.005  0.01 0.035 0.1";		# vector of displacement-cycle peaks, in terms of storey drift ratio
set Dincr [expr 0.001*$LCol];			# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis
set Fact $LCol;				# scale drift ratio by storey height for displacement cycles
set CycleType Full;				# you can do Full / Push / Half cycles with the proc
set Ncycles 1;				# specify the number of cycles at each peak

# create load pattern for lateral pushover load
set Hload [expr $Weight/2];			# define the lateral load as a proportion of the weight so that the pseudo time equals the lateral-load coefficient when using linear load pattern
set iPushNode "3 4";			# define nodes where lateral load is applied in static lateral analysis
pattern Plain 200 Linear {;			# define load pattern -- generalized
	foreach PushNode $iPushNode {
		load $PushNode $Hload 0.0 0.0 0.0 0.0 0.0
	}
}

# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

#  ---------------------------------    perform Static Cyclic Displacements Analysis
source LibGeneratePeaks.tcl
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
