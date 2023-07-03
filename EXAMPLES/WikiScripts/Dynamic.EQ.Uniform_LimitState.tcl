# --------------------------------------------------------------------------------------------------
# Example 5. 2D Frame --  Dynamic Earthquake Analysis
#                                 Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# source in procedures
source ReadSMDFile.tcl;		# procedure for reading GM file and converting it to proper format
set PI 3.141593
# Uniform Earthquake ground motion (uniform acceleration input at all support nodes)
set GMdirection 1;				# ground-motion direction
set GMfile "TCU047-N" ;			# ground-motion filenames
set GMfact 2;				# ground-motion scaling factor-1.4
#set GMdir "../GMfiles/";  # ground-motion file directory
set g 386.06
# display deformed shape:
set ViewScale 0.00004;			 # amplify display of deformed shape
#DisplayModel2D DeformedShape $ViewScale ; 	# display deformed shape, the scaling factor needs to be adjusted for each model
#recorder plot $dataDir/Drift.out Drift 550 10 700 410 -columns 1 2; # a window to plot the nodal displacements versus time
#recorder plot $dataDir/RX.out RX 10 400 345 400 -columns 1 2; # a window to plot the nodal displacements versus time

# set up ground-motion-analysis parameters
set DtAnalysis	[expr 0.0008];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr  11.5];	# maximum duration of ground-motion analysis -- should be 50*$sec

# ----------- set up analysis parameters
source LibAnalysisDynamicParameters_LimitState.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

# ------------ define & apply damping
# RAYLEIGH damping parameters, Where to put M/K-prop damping, switches (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/1099.htm)
#          D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
set xDamp 0.02;					# damping ratio
set MpropSwitch 1.0;
set KcurrSwitch 0.0;
set KcommSwitch 1.0;
set KinitSwitch 0.0;
set nEigenI 1;		# mode 1
set lambdaN [eigen [expr $nEigenI]];			# eigenvalue analysis for nEigenJ modes
set lambdaI [lindex $lambdaN [expr $nEigenI-1]]; 		# eigenvalue mode i
set omegaI [expr pow($lambdaI,0.5)];
set T1 [expr 2*$PI/$omegaI]
puts "T1"
puts "$T1 Sec"
set alphaM [expr $MpropSwitch*$xDamp*(2*$omegaI)/($omegaI)];	# M-prop. damping; D = alphaM*M
set betaKcurr [expr $KcurrSwitch*2.*$xDamp/($omegaI)];         		# current-K;      +beatKcurr*KCurrent
set betaKcomm [expr $KcommSwitch*2.*$xDamp/($omegaI)];   		# last-committed K;   +betaKcomm*KlastCommitt
set betaKinit [expr $KinitSwitch*2.*$xDamp/($omegaI)];         			# initial-K;     +beatKinit*Kini
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 				# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
set Tol 1.0e-8
# the following commands are unique to the Uniform Earthquake excitation
set IDloadTag 999;	# for uniformSupport excitation
# Uniform EXCITATION: acceleration input
set inFile $GMfile.tcl
set outFile $GMfile.g3;	# set variable holding new filename (PEER files have .at2/dt2 extension)
ReadSMDFile $inFile $outFile dt;		# call procedure to convert the ground-motion file
puts "dt is $dt"
set GMfatt [expr $g*$GMfact];		# data in input file is in g Unifts -- ACCELERATION TH
set AccelSeries "Series -dt $dt -filePath $outFile -factor  $GMfatt";	# time series information
pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $AccelSeries  ;		# create Unifform excitation

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
			test NormDispIncr   $Tol 1000  1
			algorithm ModifiedNewton -initial
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
