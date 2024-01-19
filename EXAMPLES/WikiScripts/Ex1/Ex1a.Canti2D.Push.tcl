# --------------------------------------------------------------------------------------------------
# Example 1. cantilever 2D
# static pushover analysis with gravity.
# all units are in kip, inch, second
# elasticBeamColumn ELEMENT
#			Silvia Mazzoni & Frank McKenna, 2006
#
#    ^Y
#    |
#    2       __ 
#    |         | 
#    |         | 
#    |         | 
#  (1)      36'
#    |         | 
#    |         | 
#    |         | 
#  =1=    ----  -------->X
#

# SET UP ----------------------------------------------------------------------------
wipe;						# clear opensees model
model basic -ndm 2 -ndf 3;				# 2 dimensions, 3 dof per node
file mkdir Data; 					# create data directory

# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1 0 0;					# node#, X Y
node 2 0 432

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			# node DX DY RZ

# nodal masses:
mass 2 5.18 0. 0.;					# node#, Mx My Mz, Mass=Weight/g.

# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
geomTransf Linear 1;  		# associate a tag to transformation

# connectivity: (make A very large, 10e6 times its actual value)
element elasticBeamColumn 1 1 2 3600000000 4227 1080000 1;	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag

# Define RECORDERS -------------------------------------------------------------
recorder Node -file Data/DFree.out -time -node 2 -dof 1 2 3 disp;			# displacements of free nodes
recorder Node -file Data/DBase.out -time -node 1 -dof 1 2 3 disp;			# displacements of support nodes
recorder Node -file Data/RBase.out -time -node 1 -dof 1 2 3 reaction;			# support reaction
recorder Drift -file Data/Drift.out -time -iNode 1 -jNode 2 -dof 1  -perpDirn 2 ;		# lateral drift
recorder Element -file Data/FCol.out -time -ele 1 globalForce;			# element forces -- column
recorder Element -file Data/DCol.out -time -ele 1 deformation;			# element deformations -- column

# define GRAVITY -------------------------------------------------------------
pattern Plain 1 Linear {
   load 2 0. -2000. 0.;				# node#, FX FY MZ --  superstructure-weight
}
constraints Plain;     				# how it handles boundary conditions
numberer Plain;					# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;				# how to store and solve the system of equations in the analysis
test NormDispIncr 1.0e-8 6 ; 				# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
integrator LoadControl 0.1;				# determine the next time step for an analysis, # apply gravity in 10 steps
analysis Static					# define type of analysis static or transient
analyze 10;					# perform gravity analysis
loadConst -time 0.0;				# hold gravity constant and restart time

# define LATERAL load -------------------------------------------------------------
# Lateral load pattern
pattern Plain 2 Linear {
	load 2 2000. 0.0 0.0;			# node#, FX FY MZ -- representative lateral load at top node
}

# pushover: diplacement controlled static analysis
integrator DisplacementControl 2 1 0.1;		# switch to displacement control, for node 11, dof 1, 0.1 increment
analyze 1000;				# apply 100 steps of pushover analysis to a displacement of 10

puts "Done!"




