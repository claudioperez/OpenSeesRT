# --------------------------------------------------------------------------------------------------
# Example 8. Bidirectional Uniform Eartquake Excitation
#                             Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# source in procedures
source ReadSMDFile_IR.tcl;		# procedure for reading GM file and converting it to proper format

# Bidirectional Uniform Earthquake ground motion (uniform acceleration input at all support nodes)
set iGMfile "PUL104 PUL194" ;		# ground-motion filenames, should be different files
set iGMdirection "1 3";			# ground-motion direction
set iGMfact "1.0 1.0";			# ground-motion scaling factor

# Define DISPLAY -------------------------------------------------------------
# the deformed shape is defined in the build file
#recorder plot $dataDir/DFree.out DisplDOF[lindex $iGMdirection 0] 1200 10 400 400 -columns  1 [expr 1+[lindex $iGMdirection 0]] ; # a window to plot the nodal displacements versus time
#recorder plot $dataDir/DFree.out DisplDOF[lindex $iGMdirection 1] 1200 410 400 400 -columns 1 [expr 1+[lindex $iGMdirection 1]] ; # a window to plot the nodal displacements versus time

# set up ground-motion-analysis parameters
set DtAnalysis	[expr 0.02*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr 40. *$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec


# ----------- set up analysis parameters
source LibAnalysisDynamicParameters_IR.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

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
set lambdaN [eigen -fullGenLapack [expr $nEigenJ]];			# eigenvalue analysis for nEigenJ modes
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
# the following commands are unique to the Uniform Earthquake excitation
set IDloadTag 400;	# for uniformSupport excitation
# Uniform EXCITATION: acceleration input
foreach GMdirection $iGMdirection GMfile $iGMfile GMfact $iGMfact {
	incr IDloadTag;
	set inFile $GMdir/$GMfile.tcl
	set outFile $GMdir/$GMfile.g3;			# set variable holding new filename (PEER files have .at2/dt2 extension)
	ReadSMDFile $inFile $outFile dt;			# call procedure to convert the ground-motion file
	set GMfatt [expr $g*$GMfact];			# data in input file is in g Unifts -- ACCELERATION TH
	set AccelSeries "Series -dt $dt -filePath $outFile -factor  $GMfatt";		# time series information
	pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $AccelSeries  ;	# create Unifform excitation
}

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

#