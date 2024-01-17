##########################################################
# LibMaterialsRC.tcl:  define a library of Reinforced-Concrete materials
#			Silvia Mazzoni & Frank McKenna, 2006
##########################################################

# General Material parameters
set G $Ubig;		# make stiff shear modulus
set J 1.0;			# torsional section stiffness (G makes GJ large)
set GJ [expr $G*$J];

# -----------------------------------------------------------------------------------------------------# confined and unconfined CONCRETE
# nominal concrete compressive strength
set fc 		[expr -4.0*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec 		[expr 57*$ksi*sqrt(-$fc/$psi)];	# Concrete Elastic Modulus
set nu 0.2;
set Gc [expr $Ec/2./[expr 1+$nu]];  	# Torsional stiffness Modulus

# confined concrete
set Kfc 1.3;			# ratio of confined to unconfined concrete strength
set Kres 0.2;			# ratio of residual/ultimate to maximum stress
set fc1C [expr $Kfc*$fc];		# CONFINED concrete (mander model), maximum stress
set eps1C [expr 2.*$fc1C/$Ec];	# strain at maximum stress 
set fc2C [expr $Kres*$fc1C];		# ultimate stress
set eps2C  [expr 20*$eps1C];		# strain at ultimate stress 
set lambda 0.1;			# ratio between unloading slope at $eps2 and initial slope $Ec
# unconfined concrete
set fc1U  $fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U -0.003;			# strain at maximum strength of unconfined concrete
set fc2U [expr $Kres*$fc1U];		# ultimate stress
set eps2U -0.01;			# strain at ultimate stress

# tensile-strength properties
set ftC [expr -0.14*$fc1C];		# tensile strength +tension
set ftU [expr -0.14*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness

# set up library of materials
if {  [info exists imat ] != 1} {set imat 0};		# set value only if it has not been defined previously.
set IDconcCore 1
set IDconcCover 2
uniaxialMaterial Concrete02 $IDconcCore $fc1C $eps1C $fc2C $eps2C $lambda $ftC $Ets;	# Core concrete (confined)
uniaxialMaterial Concrete02 $IDconcCover $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# Cover concrete (unconfined)

# -----------------------------------------------------------------------------------------------------# REINFORCING STEEL parameters
#
set Fy [expr 66.8*$ksi];		# STEEL yield stress
set Es [expr 29000.*$ksi];		# modulus of steel
set Bs 0.01;			# strain-hardening ratio 
set R0 18;			# control the transition from elastic to plastic branches
set cR1 0.925;			# control the transition from elastic to plastic branches
set cR2 0.15;			# control the transition from elastic to plastic branches

set IDSteel 3
uniaxialMaterial Steel02 $IDSteel  $Fy $Es $Bs $R0 $cR1 $cR2
