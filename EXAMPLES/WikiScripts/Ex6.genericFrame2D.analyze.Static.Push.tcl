# --------------------------------------------------------------------------------------------------
# Example 6. 2D Frame --  Static Pushover Analysis
#                                Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# characteristics of pushover analysis
set Dmax [expr 0.1*$LBuilding ];	# maximum displacement of pushover. push to 10% drift.
set Dincr [expr 0.001*$LBuilding ];	# displacement increment. you want this to be small, but not too small to slow analysis

# -- STATIC PUSHOVER/CYCLIC ANALYSIS
# create load pattern for lateral pushover load coefficient when using linear load pattern
pattern Plain 200 Linear {;			# define load pattern
	foreach FPush $iFPush NodePush $iNodePush {;		# these vectors have been defined with the structure
			load $NodePush  $FPush 0.0 0.0
	}
};		# end load pattern

# display deformed shape:
set ViewScale 5;
#DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model
#recorder plot $dataDir/DFree.out ForceDisp 910 10 400 400 -columns 2 1; # a window to plot the load vs. nodal displacement

#  ---------------------------------    perform Static Pushover Analysis
# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis
# ----------------------------------------------first analyze command------------------------
set Nsteps [expr int($Dmax/$Dincr)];		# number of pushover analysis steps
set ok [analyze $Nsteps];                		# this will return zero if no convergence problems were encountered
# ----------------------------------------------if convergence failure-------------------------
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
