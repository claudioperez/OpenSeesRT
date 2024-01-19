# --------------------------------------------------------------------------------------------------
# Example 3. 2D Cantilever -- EQ ground motion
#                                 Silvia Mazzoni & Frank McKenna, 2006
# execute this file after you have built the model, and after you apply gravity
#

# Uniform Earthquake ground motion (uniform acceleration input at all support nodes)
set GMdirection 1;				# ground-motion direction
set GMfile "BM68elc.acc" ;			# ground-motion filenames
set GMfact 1.5;				# ground-motion scaling factor

# set up ground-motion-analysis parameters
set DtAnalysis	[expr 0.01*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr 10. *$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec

# DYNAMIC ANALYSIS PARAMETERS
# CONSTRAINTS handler -- Determines how the constraint equations are enforced in the analysis (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/617.htm)
#          Plain Constraints -- Removes constrained degrees of freedom from the system of equations 
#          Lagrange Multipliers -- Uses the method of Lagrange multipliers to enforce constraints 
#          Penalty Method -- Uses penalty numbers to enforce constraints 
#          Transformation Method -- Performs a condensation of constrained degrees of freedom 
constraints Transformation ; 

# DOF NUMBERER (number the degrees of freedom in the domain): (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/366.htm)
#   determines the mapping between equation numbers and degrees-of-freedom
#          Plain -- Uses the numbering provided by the user 
#          RCM -- Renumbers the DOF to minimize the matrix band-width using the Reverse Cuthill-McKee algorithm 
numberer Plain

# SYSTEM (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/371.htm)
#   Linear Equation Solvers (how to store and solve the system of equations in the analysis)
#   -- provide the solution of the linear system of equations Ku = P. Each solver is tailored to a specific matrix topology. 
#          ProfileSPD -- Direct profile solver for symmetric positive definite matrices 
#          BandGeneral -- Direct solver for banded unsymmetric matrices 
#          BandSPD -- Direct solver for banded symmetric positive definite matrices 
#          SparseGeneral -- Direct solver for unsymmetric sparse matrices (-piv option)
#          SparseSPD -- Direct solver for symmetric sparse matrices 
#          UmfPack -- Direct UmfPack solver for unsymmetric matrices 
system SparseGeneral -piv

# TEST: # convergence test to 
# Convergence TEST (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/360.htm)
#   -- Accept the current state of the domain as being on the converged solution path 
#   -- determine if convergence has been achieved at the end of an iteration step
#          NormUnbalance -- Specifies a tolerance on the norm of the unbalanced load at the current iteration 
#          NormDispIncr -- Specifies a tolerance on the norm of the displacement increments at the current iteration 
#          EnergyIncr-- Specifies a tolerance on the inner product of the unbalanced load and displacement increments at the current iteration 
#          RelativeNormUnbalance --
#          RelativeNormDispIncr --
#          RelativeEnergyIncr --
set Tol 1.e-8;                        # Convergence Test: tolerance
set maxNumIter 10;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlag 0;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set TestType EnergyIncr;	# Convergence-test type
test $TestType $Tol $maxNumIter $printFlag;

# Solution ALGORITHM: -- Iterate from the last time step to the current (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/682.htm)
#          Linear -- Uses the solution at the first iteration and continues 
#          Newton -- Uses the tangent at the current iteration to iterate to convergence 
#          ModifiedNewton -- Uses the tangent at the first iteration to iterate to convergence 
#          NewtonLineSearch -- 
#          KrylovNewton -- 
#          BFGS -- 
#          Broyden -- 
set algorithmType ModifiedNewton 
algorithm $algorithmType;        

# Static INTEGRATOR: -- determine the next time step for an analysis  (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/689.htm)
#          LoadControl -- Specifies the incremental load factor to be applied to the loads in the domain 
#          DisplacementControl -- Specifies the incremental displacement at a specified DOF in the domain 
#          Minimum Unbalanced Displacement Norm -- Specifies the incremental load factor such that the residual displacement norm in minimized 
#          Arc Length -- Specifies the incremental arc-length of the load-displacement path 
# Transient INTEGRATOR: -- determine the next time step for an analysis including inertial effects 
#          Newmark -- The two parameter time-stepping method developed by Newmark 
#          HHT -- The three parameter Hilbert-Hughes-Taylor time-stepping method 
#          Central Difference -- Approximates velocity and acceleration by centered finite differences of displacement 
set NewmarkGamma 0.5;	# Newmark-integrator gamma parameter (also HHT)
set NewmarkBeta 0.25;	# Newmark-integrator beta parameter
integrator Newmark $NewmarkGamma $NewmarkBeta 

# ANALYSIS  -- defines what type of analysis is to be performed (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/324.htm)
#          Static Analysis -- solves the KU=R problem, without the mass or damping matrices. 
#          Transient Analysis -- solves the time-dependent analysis. The time step in this type of analysis is constant. The time step in the output is also constant. 
#          variableTransient Analysis -- performs the same analysis type as the Transient Analysis object. The time step, however, is variable. This method is used when 
#                 there are convergence problems with the Transient Analysis object at a peak or when the time step is too small. The time step in the output is also variable.
analysis Transient

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
# define damping
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 				# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# Uniform EXCITATION: acceleration input
set IDloadTag 400;			# load tag
set dt 0.01;			# time step for input ground motion
set GMfatt 1.0;			# data in input file is in g Unifts -- ACCELERATION TH
set AccelSeries "Series -dt $dt -filePath $GMfile -factor  $GMfatt";			# time series information
pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $AccelSeries  ;		# create Unifform excitation

set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];
set ok [analyze $Nsteps $DtAnalysis];			# actually perform analysis; returns ok=0 if analysis was successful

if {$ok != 0} {      ;					# if analysis was not successful.
	# change some analysis parameters to achieve convergence
	# performance is slower inside this loop
	#    Time-controlled analysis
	set ok 0;
	set controlTime [getTime];
	while {$controlTime < $TmaxAnalysis && $ok == 0} {
		set ok [analyze 1 $DtAnalysis]
		set controlTime [getTime]
		set ok [analyze 1 $DtAnalysis]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 1000  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $TestType $Tol $maxNumIter  0
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmType
		}
	}
};      # end if ok !0



puts "Ground Motion Done. End Time: [getTime]"
