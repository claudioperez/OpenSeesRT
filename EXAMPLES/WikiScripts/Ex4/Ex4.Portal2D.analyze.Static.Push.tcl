# --------------------------------------------------------------------------------------------------
# Example4. 2D Portal Frame--  Static Pushover Analysis
#                            Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# we need to set up parameters that are particular to the model.
set IDctrlNode 3;			# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement contro
# characteristics of pushover analysis
set Dmax [expr 0.1*$LCol];		# maximum displacement of pushover. push to 10% drift.
set Dincr [expr 0.001*$LCol];		# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis

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

#  ---------------------------------    perform Static Pushover Analysis
set Nsteps [expr int($Dmax/$Dincr)];        # number of pushover analysis steps
set ok [analyze $Nsteps];                # this will return zero if no convergence problems were encountered
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis
if {$ok != 0} {  
	# if analysis fails, we try some other stuff, performance is slower inside this loop
	set Dstep 0.0;
	set ok 0
	while {$Dstep <= 1.0 && $ok == 0} {	
		set controlDisp [nodeDisp $IDctrlNode $IDctrlDOF ]
		set Dstep [expr $controlDisp/$Dmax]
		set ok [analyze 1 ]
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

	};	# end while loop
};      # end if ok !0

# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}
