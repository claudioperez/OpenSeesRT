# --------------------------------------------------------------------------------------------------
# Example 7. Static Pushover
#                               Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# characteristics of pushover analysis
set Dmax [expr 0.1*$LBuilding ];	# maximum displacement of pushover. push to 10% drift.
set Dincr [expr 0.001*$LBuilding ];	# displacement increment. you want this to be small, but not too small to slow analysis

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
set ViewScale 1;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels
#recorder plot $dataDir/DFree.out Displ-X [expr $xPixels+10] 10 300 300 -columns 2 1; # a window to plot the nodal displacements versus time

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
		set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
		if {$ok != 0} {;				# reduce step size if still fails to converge
			set Nk 4;			# reduce step size
			set DincrReduced [expr $Dincr/$Nk];
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $DincrReduced
			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
				set ok [analyze 1];                		# this will return zero if no convergence problems were encountered
				if {$ok != 0} {  
					# if analysis fails, we try some other stuff
					# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
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
				if {$ok != 0} {;				# stop if still fails to converge
					puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					return -1
				}; # end if
			}; # end for
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;	# bring back to original increment
		}; # end if
	};	# end while loop
};      # end if ok !0
# -----------------------------------------------------------------------------------------------------

if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}
