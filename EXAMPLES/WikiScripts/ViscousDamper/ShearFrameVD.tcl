# #####################################################################################
#
# Modelling of Single Story Shear Frame equipped with Nonlinear Viscous Damper
#
# Sarven Akcelyan & Dimitrios G. Lignos
# McGill University, Quebec, Canada
#
# Date: 20/08/2013
# Revised: 20/08/2013
#
# #####################################################################################
# Define model
# All the units are in mm,KN,sec

wipe all;							# clear memory of past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm = #dimension, ndf = #dofs

# create data directory
set Output Output;
file mkdir $Output; 

# define geometry
set L  5000.; #bay width
set h  3000.; #story heigth

# define nodal coordinates:
node 1   0.  0. ;
node 2   $L  0. ;
node 3   0.  $h ;
node 4   $L  $h ;

# Single point constraints -- Boundary Conditions
fix 1 1 1 1;
fix 2 1 1 1;

# MP constraints
equalDOF 3 4 2 3 ; #Shear Beam

# mass
set W 1000.; #KN
set g 9810.; #mm/sec2
set m  [expr $W/$g]; 

# assign mass
mass 3 [expr 0.5*$m] 0. 0. ;
mass 4 [expr 0.5*$m] 0. 0. ;
				
				
set Tn 0.7; # sec (Natural Period)
set pi [expr acos(-1.0)];

# Columns and Beam Properties
set K [expr pow(2*$pi/$Tn,2)*$m]; # KN/mm
set E 200.0; # KN/mm2
set Ic [expr $K*pow($h,3)/(24*$E)]; # mm4 (K=24EIc/h^3)
set Ib [expr 1e12*$Ic];  # mm4
set A  [expr 1e12]; # mm2.

# Damper Properties
set Kd 25.
set Cd 20.7452;
set ad 0.35;

# Define ViscousDamper Material
#uniaxialMaterial ViscousDamper $matTag $Kd $Cd $alpha
uniaxialMaterial ViscousDamper     1   $Kd  $Cd $ad

# define geometric transformation:
set TransfTag 1;
geomTransf Linear $TransfTag ;

# Define Elements------------------------------------------------------------
# Columns
element elasticBeamColumn 1 1 3 $A $E $Ic $TransfTag;
element elasticBeamColumn 2 2 4 $A $E $Ic $TransfTag;
# Beam
element elasticBeamColumn 3 3 4 $A $E $Ib $TransfTag;
# Damper
#element twoNodeLink $eleTag $iNode $jNode -mat $matTags -dir $dirs
element  twoNodeLink   4      1      4     -mat   1      -dir 1 

puts "Model Built"
	
# Record--------------------------------------------------------------------

#timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor>
timeSeries Path 1 -dt 0.01 -filePath TakY.th -factor [expr 0.50*$g];  # define acceleration vector from file (dt=0.01 is associated with the input file gm)

# Define RECORDERS -------------------------------------------------------------
recorder Node -file $Output/Disp.out -time -node 4 -dof 1  disp;
recorder Node -file $Output/Acc.out -timeSeries 1 -time -node 4 -dof 1  accel;

recorder Node -file $Output/Base.out -time -node 1 2 -dof 1  reaction;  # support reaction
recorder Node -file $Output/NBase.out -time -node 1 2 -dof 2  reaction;	 # support reaction
			       
recorder Element -file  $Output/Damperdisp.out -time  -ele 4 deformations ; 
recorder Element -file  $Output/Damperforce.out -time -ele 4 localForce ;
recorder Element -file  $Output/Dampergbforce.out -time -ele 4 -dof 1  force;
recorder Element -file  $Output/Frameforce.out -time -ele 1 2 -dof 1  force;

# set damping based on first eigen mode
set freq [expr [eigen 1]**0.5]
set period [expr 2*$pi/$freq]
puts $period
set damp 0.02;
rayleigh [expr 2*$damp*$freq] 0.  0. 0.   

#pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $ver0>
pattern UniformExcitation 1 1 -accel 1;		         # define where and how (pattern tag, dof) acceleration is applied


# display displacement shape of the column
# recorder display "Displaced shape" 10 10 500 500 -wipe
# prp 200. 50. 1;
# vup  0  1 0;
# vpn  0  0 1;
# display 1 5 40 

# create the analysis
wipeAnalysis;					     # clear previously-define analysis parameters
constraints Transformation;     	 # how it handles boundary conditions
numberer RCM;					     # renumber dof's to minimize band-width (optimization), if you want to
system UmfPack;						 # how to store and solve the system of equations in the analysis (large model: try UmfPack)
test EnergyIncr 1.0e-10 100; 		 # test Eneregy incerment 
algorithm  KrylovNewton;		 	 # use Kyrlow-Newton algorithm 
integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
analysis Transient;					 # define type of analysis: time-dependent
analyze [expr 10*4096] 0.001;		 # apply 10*4096 steps for 0.001-sec time steps in analysis

puts "Done!"
wipe
