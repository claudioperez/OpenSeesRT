# --------------------------------------------------------------------------------------------------
# dynamic-analysis parameters
# I am setting all these variables as global variables (using variable rather than set command)
#    so that these variables can be uploaded by a procedure
#                                 Silvia Mazzoni & Frank McKenna, 2006


# Set up Analysis Parameters ---------------------------------------------
# CONSTRAINTS handler -- Determines how the constraint equations are enforced in the analysis (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/617.htm)
#          Plain Constraints -- Removes constrained degrees of freedom from the system of equations 
#          Lagrange Multipliers -- Uses the method of Lagrange multipliers to enforce constraints 
#          Penalty Method -- Uses penalty numbers to enforce constraints 
#          Transformation Method -- Performs a condensation of constrained degrees of freedom 
variable constraintsTypeDynamic Transformation;
constraints $constraintsTypeDynamic ; 

# DOF NUMBERER (number the degrees of freedom in the domain): (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/366.htm)
#   determines the mapping between equation numbers and degrees-of-freedom
#          Plain -- Uses the numbering provided by the user 
#          RCM -- Renumbers the DOF to minimize the matrix band-width using the Reverse Cuthill-McKee algorithm 
variable numbererTypeDynamic RCM
numberer $numbererTypeDynamic 

# SYSTEM (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/371.htm)
#   Linear Equation Solvers (how to store and solve the system of equations in the analysis)
#   -- provide the solution of the linear system of equations Ku = P. Each solver is tailored to a specific matrix topology. 
#          ProfileSPD -- Direct profile solver for symmetric positive definite matrices 
#          BandGeneral -- Direct solver for banded unsymmetric matrices 
#          BandSPD -- Direct solver for banded symmetric positive definite matrices 
#          SparseGeneral -- Direct solver for unsymmetric sparse matrices (-piv option)
#          SparseSPD -- Direct solver for symmetric sparse matrices 
#          UmfPack -- Direct UmfPack solver for unsymmetric matrices 
variable systemTypeDynamic BandGeneral;	# try UmfPack for large problems
system $systemTypeDynamic 

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
variable TolDynamic 1.e-8;                        # Convergence Test: tolerance
variable maxNumIterDynamic 10;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
variable printFlagDynamic 0;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
variable testTypeDynamic EnergyIncr;	# Convergence-test type
test $testTypeDynamic $TolDynamic $maxNumIterDynamic $printFlagDynamic;
# for improved-convergence procedure:
	variable maxNumIterConvergeDynamic 2000;	
	variable printFlagConvergeDynamic 0;	

# Solution ALGORITHM: -- Iterate from the last time step to the current (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/682.htm)
#          Linear -- Uses the solution at the first iteration and continues 
#          Newton -- Uses the tangent at the current iteration to iterate to convergence 
#          ModifiedNewton -- Uses the tangent at the first iteration to iterate to convergence 
#          NewtonLineSearch -- 
#          KrylovNewton -- 
#          BFGS -- 
#          Broyden -- 
variable algorithmTypeDynamic ModifiedNewton 
algorithm $algorithmTypeDynamic;        

# Static INTEGRATOR: -- determine the next time step for an analysis  (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/689.htm)
#          LoadControl -- Specifies the incremental load factor to be applied to the loads in the domain 
#          DisplacementControl -- Specifies the incremental displacement at a specified DOF in the domain 
#          Minimum Unbalanced Displacement Norm -- Specifies the incremental load factor such that the residual displacement norm in minimized 
#          Arc Length -- Specifies the incremental arc-length of the load-displacement path 
# Transient INTEGRATOR: -- determine the next time step for an analysis including inertial effects 
#          Newmark -- The two parameter time-stepping method developed by Newmark 
#          HHT -- The three parameter Hilbert-Hughes-Taylor time-stepping method 
#          Central Difference -- Approximates velocity and acceleration by centered finite differences of displacement 
variable NewmarkGamma 0.5;	# Newmark-integrator gamma parameter (also HHT)
variable NewmarkBeta 0.25;	# Newmark-integrator beta parameter
variable integratorTypeDynamic Newmark;
integrator $integratorTypeDynamic $NewmarkGamma $NewmarkBeta

# ANALYSIS  -- defines what type of analysis is to be performed (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/324.htm)
#          Static Analysis -- solves the KU=R problem, without the mass or damping matrices. 
#          Transient Analysis -- solves the time-dependent analysis. The time step in this type of analysis is constant. The time step in the output is also constant. 
#          variableTransient Analysis -- performs the same analysis type as the Transient Analysis object. The time step, however, is variable. This method is used when 
#                 there are convergence problems with the Transient Analysis object at a peak or when the time step is too small. The time step in the output is also variable.
variable analysisTypeDynamic Transient
analysis $analysisTypeDynamic 
