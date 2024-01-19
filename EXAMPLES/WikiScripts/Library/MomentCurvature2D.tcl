proc MomentCurvature2D { secTag axialLoad maxK {numIncr 100} } {
	##################################################
	# A procedure for performing section analysis (only does
	# moment-curvature, but can be easily modified to do any mode
	# of section reponse.)
	#
	# MHS
	# October 2000
	# modified to 2D and to improve convergence by Silvia Mazzoni, 2006
	#
	# Arguments
	#	secTag -- tag identifying section to be analyzed
	#	axialLoad -- axial load applied to section (negative is compression)
	#	maxK -- maximum curvature reached during analysis
	#	numIncr -- number of increments used to reach maxK (default 100)
	#
	# Sets up a recorder which writes moment-curvature results to file
	# section$secTag.out ... the moment is in column 1, and curvature in column 2

	# Define two nodes at (0,0)
	node 1001 0.0 0.0
	node 1002 0.0 0.0

	# Fix all degrees of freedom except axial and bending
	fix 1001 1 1 1
	fix 1002 0 1 0

	# Define element
	#                         tag ndI ndJ  secTag
	element zeroLengthSection  2001   1001   1002  $secTag

	# Create recorder
	recorder Node -file data/Mphi.out -time -node 1002 -dof 3 disp;	# output moment (col 1) & curvature (col 2)
	
	# Define constant axial load
	pattern Plain 3001 "Constant" {
		load 1002 $axialLoad 0.0 0.0
	}

	# Define analysis parameters
	integrator LoadControl 0 1 0 0
	system SparseGeneral -piv;	# Overkill, but may need the pivoting!
	test EnergyIncr  1.0e-9 10
	numberer Plain
	constraints Plain
	algorithm Newton
	analysis Static

	# Do one analysis for constant axial load
	analyze 1
	
	loadConst -time 0.0

	# Define reference moment
	pattern Plain 3002 "Linear" {
		load 1002 0.0 0.0 1.0
	}

	# Compute curvature increment
	set dK [expr $maxK/$numIncr]

	# Use displacement control at node 1002 for section analysis, dof 3
	integrator DisplacementControl 1002 3 $dK 1 $dK $dK

	# Do the section analysis
	set ok [analyze $numIncr]

	# ----------------------------------------------if convergence failure-------------------------
	set IDctrlNode 1002
	set IDctrlDOF 3
	set Dmax $maxK
	set Dincr $dK
	set TolStatic 1.e-9;
	set testTypeStatic EnergyIncr  
	set maxNumIterStatic 6
	set algorithmTypeStatic Newton
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
						test NormDispIncr   $TolStatic      2000 0
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
	global LunitTXT;					# load time-unit text
	if {  [info exists LunitTXT] != 1} {set LunitTXT "Length"};		# set blank if it has not been defined previously.

	set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Curv=%.4f /%s";	# format for screen/file output of DONE/PROBLEM analysis
	if {$ok != 0 } {
		puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
	} else {
		puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
	}

}
