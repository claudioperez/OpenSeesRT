# --------------------------------------------------------------------------------------------------
# Example 5. 2D Frame --  Dynamic Sine-Wave Input Analysis, multiple-support excitation
#                                 Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# MultipleSupport SineWave ground motion (different displacement input at spec'd support nodes) -- four nodes here
# Sine Input:
set iSupportNode "11 12 13 14" ;			# support nodes where ground motions are input, for multiple-support excitation
set iGMdirection "1 1 1 1";		# ground-motion direction  -- for each support node
set iGMSineDispAmpl "[expr 0.01*$in] [expr 0.1*$in] [expr 0.2*$in] [expr 0.3*$in]"; # sine ground-motion displacement amplitude (this is the support motion, not the free-node motion)
set iTPeriodSine "[expr 0.35*$sec] [expr 0.36*$sec] [expr 0.37*$sec] [expr 0.38*$sec] "; 	# period of input sine wave
set iDurationSine "[expr 3.*$sec] [expr 3.1*$sec] [expr 3.2*$sec] [expr 3.3*$sec] ";		# duration of input sine wave

# set up ground-motion-analysis parameters
set DtAnalysis [expr 0.01*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis [expr 10. *$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec

# ----------- set up analysis parameters
source LibAnalysisDynamicParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

# ------------ define & apply damping
# RAYLEIGH damping parameters, Where to put M/K-prop damping, switches (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/1099.htm)
#          D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
set xDamp 0.02;					# damping ratio
set MpropSwitch 1.0;
set KcurrSwitch 0.0;
set KcommSwitch 1.0;
set KinitSwitch 0.0;
set nEigenI 1;		# mode 1
set nEigenJ 3;		# mode 3
set lambdaN [eigen [expr $nEigenJ]];			# eigenvalue analysis for nEigenJ modes
set lambdaI [lindex $lambdaN [expr $nEigenI-1]]; 		# eigenvalue mode i
set lambdaJ [lindex $lambdaN [expr $nEigenJ-1]]; 	# eigenvalue mode j
set omegaI [expr pow($lambdaI,0.5)];
set omegaJ [expr pow($lambdaJ,0.5)];
set alphaM [expr $MpropSwitch*$xDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	# M-prop. damping; D = alphaM*M
set betaKcurr [expr $KcurrSwitch*2.*$xDamp/($omegaI+$omegaJ)];         		# current-K;      +beatKcurr*KCurrent
set betaKcomm [expr $KcommSwitch*2.*$xDamp/($omegaI+$omegaJ)];   		# last-committed K;   +betaKcomm*KlastCommitt
set betaKinit [expr $KinitSwitch*2.*$xDamp/($omegaI+$omegaJ)];         			# initial-K;     +beatKinit*Kini
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 				# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# the following commands are unique to the Multiple-Support Earthquake excitation
set IDloadTag 400;	
set IDgmSeries 500;	# for multipleSupport Excitation
set DtGround [expr 0.02*$sec];	# time-step Dt for input grond motion
# multiple-support excitation: displacement input at individual nodes
pattern MultipleSupport $IDloadTag  {
	foreach SupportNode $iSupportNode GMdirection $iGMdirection GMSineDispAmpl $iGMSineDispAmpl TPeriodSine $iTPeriodSine DurationSine $iDurationSine {
		set IDgmSeries [expr $IDgmSeries +1]
		set DispSeries "Sine 0. $DurationSine $TPeriodSine -factor $GMSineDispAmpl"
		groundMotion $IDgmSeries Plain -disp  $DispSeries  
	     	imposedMotion $SupportNode  $GMdirection $IDgmSeries
	};	# end foreach	
};	# end pattern

set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];
set ok [analyze $Nsteps $DtAnalysis];			# actually perform analysis; returns ok=0 if analysis was successful

if {$ok != 0} {      ;					# analysis was not successful.
	# --------------------------------------------------------------------------------------------------
	# change some analysis parameters to achieve convergence
	# performance is slower inside this loop
	#    Time-controlled analysis
	set ok 0;
	set controlTime [getTime];
	while {$controlTime < $TmaxAnalysis && $ok == 0} {
		set controlTime [getTime]
		set ok [analyze 1 $DtAnalysis]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 1000  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmTypeDynamic
		}
	}
};      # end if ok !0


puts "Ground Motion Done. End Time: [getTime]"
