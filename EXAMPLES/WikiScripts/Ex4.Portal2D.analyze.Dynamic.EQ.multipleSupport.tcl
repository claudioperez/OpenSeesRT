# --------------------------------------------------------------------------------------------------
# Example4. 2D Portal Frame--  Dynamic EQ input analysis -- multiple-support excitation
#                                Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# MultipleSupport Earthquake ground motion (different displacement input at spec'd support nodes) -- two nodes here
set iSupportNode "1 2";			# support nodes where ground motions are input, for multiple-support excitation
set iGMfact "1.5 1.25";			# ground-motion scaling factor (units not included here)  -- for each support node
set iGMdirection "1 1";			# ground-motion direction  -- for each support node
# ground-motion input:
set iGMfile "H-e12140  H-e12140";		# groun-motion filename for each support node  -- each should be different

# set up ground-motion-analysis parameters
set DtAnalysis	[expr 0.01*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr 10. *$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec

# ----------- set up analysis parameters
source LibAnalysisDynamicParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

# define DAMPING--------------------------------------------------------------------------------------
# apply Rayleigh DAMPING from $xDamp
# D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
set xDamp 0.02;				# 2% damping ratio
set lambda [eigen 1]; 			# eigenvalue mode 1
set omega [expr pow($lambda,0.5)];
set alphaM 0.;				# M-prop. damping; D = alphaM*M
set betaKcurr 0.;         			# K-proportional damping;      +beatKcurr*KCurrent
set betaKcomm [expr 2.*$xDamp/($omega)];   	# K-prop. damping parameter;   +betaKcomm*KlastCommitt
set betaKinit 0.;         			# initial-stiffness proportional damping      +beatKinit*Kini
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 		# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# the following commands are unique to the Multiple-Support Earthquake excitation
set IDloadTag 400;	
set IDgmSeries 500;	# for multipleSupport Excitation
# read a PEER strong motion database file, extracts dt from the header and converts the file 
# to the format OpenSees expects for Uniform/multiple-support ground motions 
source ReadSMDFile.tcl;	# read in procedure Multinition
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