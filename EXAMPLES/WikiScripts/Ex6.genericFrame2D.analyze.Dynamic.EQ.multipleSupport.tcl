# --------------------------------------------------------------------------------------------------
# Example 6. 2D Frame --  Multiple-Support EQ ground motion
#                                 Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# source in procedures
source ReadSMDFile.tcl;		# procedure for reading GM file and converting it to proper format

# MultipleSupport Earthquake ground motion (different displacement input at spec'd support nodes) 
# iSupportNode, supportnodes where ground motions are applied, is defined with model
set iGMfact "1.0 1.1 1.2 1.3"
set iGMdirection "1 1 1 1"
set iGMfile "H-e12140 H-e12140 H-e12140 H-e12140"

# display deformed shape:
set ViewScale 5;
#DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model
#recorder plot $dataDir/DFree.out Displ 10 700 400 400 -columns 1 2; # a window to plot the nodal displacements versus time

# set up ground-motion-analysis parameters
set DtAnalysis	[expr 0.01*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr 10. *$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec

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

# define damping
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 				# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# the following commands are unique to the Multiple-Support Earthquake excitation
set IDloadTag 400;	
set IDgmSeries 500;	# for multipleSupport Excitation
# multiple-support excitation: displacement input at individual nodes
pattern MultipleSupport $IDloadTag  {
	foreach SupportNode $iSupportNode GMfile $iGMfile GMfact $iGMfact GMdirection $iGMdirection {
		set IDgmSeries [expr $IDgmSeries +1]
		set inFile $GMdir/$GMfile.dt2
		set outFile $GMdir/$GMfile.g3;	# set variable holding new filename (PEER files have .at2/dt2 extension)
		ReadSMDFile $inFile $outFile dt;		# call procedure to convert the ground-motion file
		set GMfatt [expr $cm*$GMfact];	# data in input file is in cm Unifts -- DISPLACEMENT TH
		set DispSeries "Series -dt $dt -filePath $outFile -factor  $GMfatt";	# time series information
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
