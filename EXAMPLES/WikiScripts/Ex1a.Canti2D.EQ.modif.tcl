# --------------------------------------------------------------------------------------------------
# Example 1. cantilever 2D
# EQ ground motion with gravity
# all units are in kip, inch, second
# elasticBeamColumn ELEMENT
#		Silvia Mazzoni & Frank McKenna, 2006
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
wipe;						       # clear opensees model
model basic -ndm 2 -ndf 3;	       # 2 dimensions, 3 dof per node
file mkdir data; 				   # create data directory

# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1 0. 0.;					   # node#, X Y
node 2 0. 432.

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			           # node DX DY RZ

# nodal masses:
mass 2 5.18 0. 0.;			   # node#, Mx My Mz, Mass=Weight/g.

# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
geomTransf Linear 1;  		       # associate a tag to transformation

# connectivity:
element elasticBeamColumn 1 1 2 3600 3225 1080000 1;	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag

# Define RECORDERS -------------------------------------------------------------
recorder Node -file Data/DFree.out -time -node 2 -dof 1 2 3 disp;			            # displacements of free nodes
recorder Node -file Data/RBase.out -time -node 1 -dof 1 2 3 reaction;			        # support reaction
recorder Drift -file Data/Drift.out -time -iNode 1 -jNode 2 -dof 1  -perpDirn 2 ;		# lateral drift
recorder Element -file Data/FCol.out -time -ele 1 force;			                    # element forces -- column

# define GRAVITY -------------------------------------------------------------
timeSeries Linear 1
pattern Plain 1 1 {
   load 2 0. -2000. 0.;			    # node#, FX FY MZ --  superstructure-weight
}
constraints Plain;     				# how it handles boundary conditions
numberer Plain;					    # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;				    # how to store and solve the system of equations in the analysis
algorithm Linear;                   # use Linear algorithm for linear analysis
integrator LoadControl 0.1;			# determine the next time step for an analysis, # apply gravity in 10 steps
analysis Static					    # define type of analysis static or transient
analyze 10;					        # perform gravity analysis
loadConst -time 0.0;				# hold gravity constant and restart time

# DYNAMIC ground-motion analysis -------------------------------------------------------------
# create load pattern
set G 386
timeSeries Path 2 -dt 0.005 -filePath A10000.tcl -factor $G; # define acceleration vector from file (dt=0.005 is associated with the input file gm)
pattern UniformExcitation 2 1 -accel 2;		         # define where and how (pattern tag, dof) acceleration is applied

# set damping based on first eigen mode
set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio 0.02
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]

# display displacement shape of the column
# recorder display "Displaced shape" 10 10 500 500 -wipe
prp 200. 50. 1;
vup  0  1 0;
vpn  0  0 1;
display 1 5 40 

# create the analysis
wipeAnalysis;					     # clear previously-define analysis parameters
constraints Plain;     				 # how it handles boundary conditions
numberer Plain;					     # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;					 # how to store and solve the system of equations in the analysis
algorithm Linear					 # use Linear algorithm for linear analysis
integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
analysis Transient;					 # define type of analysis: time-dependent
analyze 3995 0.01;					 # apply 3995 0.01-sec time steps in analysis


puts "Done!"
wipe

