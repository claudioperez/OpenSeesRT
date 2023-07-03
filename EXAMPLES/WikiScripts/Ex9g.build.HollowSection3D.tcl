# --------------------------------------------------------------------------------------------------
# build a hollow reinforced concrete confined section
#		Vesna Terzic, 2010
#

# SET UP ----------------------------------------------------------------------------
wipe;				# clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name of data directory -- simple
file mkdir $dataDir; 			# create data directory
source LibUnits.tcl;			# define units

# MATERIAL parameters -------------------------------------------------------------------
set IDconcCore 1; 				# material ID tag -- confined core concrete
set IDconcCover 2; 				# material ID tag -- unconfined cover concrete
set IDreinf 3; 				# material ID tag -- reinforcement
# nominal concrete compressive strength
set fc 		[expr -4.0*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec 		[expr 57*$ksi*sqrt(-$fc/$psi)];	# Concrete Elastic Modulus
# confined concrete
set Kfc 		1.3;			# ratio of confined to unconfined concrete strength
set fc1C 		[expr $Kfc*$fc];		# CONFINED concrete (mander model), maximum stress
set eps1C	[expr 2.*$fc1C/$Ec];	# strain at maximum stress 
set fc2C 		[expr 0.2*$fc1C];		# ultimate stress
set eps2C 	[expr 5*$eps1C];		# strain at ultimate stress 
# unconfined concrete
set fc1U 		$fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U	-0.003;			# strain at maximum strength of unconfined concrete
set fc2U 		[expr 0.2*$fc1U];		# ultimate stress
set eps2U	-0.01;			# strain at ultimate stress
set lambda 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
set ftC [expr -0.07*$fc1C];		# tensile strength +tension
set ftU [expr -0.07*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness
# -----------
set Fy 		[expr 66.8*$ksi];		# STEEL yield stress
set Es		[expr 29000.*$ksi];		# modulus of steel
set Bs		0.01;			# strain-hardening ratio 
set R0 18;				# control the transition from elastic to plastic branches
set cR1 0.925;				# control the transition from elastic to plastic branches
set cR2 0.15;				# control the transition from elastic to plastic branches
uniaxialMaterial Concrete02 $IDconcCore $fc1C $eps1C $fc2C $eps2C $lambda $ftC $Ets;	# build core concrete (confined)
uniaxialMaterial Concrete02 $IDconcCover $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# build cover concrete (unconfined)
uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;				# build reinforcement material

# Geometry of the section and reinforcement -------------------------------------------------------------------
set SecTag 1; # section ID tag
set pi 3.14
set b [expr 400*$cm]; # width of the section
set d [expr 560*$cm]; # depth of the section
set bh [expr 300*$cm]; # width of the hole
set dh [expr 460*$cm]; # depth of the hole
set cover [expr 5*$cm]; # concrete cover
set barD [expr 3.2*$cm]; # longitudianl bar diameter
set areaFiber [expr $barD**2*$pi/4]; # area of the longitudianl bar
set dStirrup [expr 1*$cm]; #diameter of stirrups 
set shift [expr 15*$cm]; #distance from the beginning of the web to the first reinf. bar

set numFiber1 25; #number of reinforceng bars in layer type 1
set numFiber2 13; #number of reinforceng bars in layer type 2

#number of subdivision for one patch 
set numSubdivIJ1 112
set numSubdivJK1 1
set numSubdivIJ2 92
set numSubdivJK2 1
set numSubdivIJ3 1
set numSubdivJK3 78
set numSubdivIJ4 1
set numSubdivJK4 62
set numSubdivIJ5 110
set numSubdivJK5 8
set numSubdivIJ6 8
set numSubdivJK6 62

#coordinates that define different patches of confined and unconfined concrete and layers of reinforcement
set y1 [expr $d/2]
set z1 [expr $b/2]
set z2 [expr $z1-$cover-$dStirrup-$barD/2]
set y2 [expr $dh/2]
set z4 [expr $bh/2]
set z3 [expr $z4+$cover+$dStirrup+$barD/2]
set y3 [expr $y1-$cover-$dStirrup-$barD/2]
set y4 [expr $y2+$cover+$dStirrup+$barD/2]

#coordiantes for steel layers 
set z4s [expr $z4-$shift]

section Fiber $SecTag {
	#cover
	#                          nfIJ              nfJK    yI  zI  yJ  zJ  yK  zK    yL  zL
	patch quad $IDconcCover $numSubdivIJ1 $numSubdivJK1 -$y1 $z2 $y1 $z2 $y1 $z1 -$y1 $z1
	patch quad $IDconcCover $numSubdivIJ1 $numSubdivJK1 -$y1 -$z1 $y1 -$z1 $y1 -$z2 -$y1 -$z2
	patch quad $IDconcCover $numSubdivIJ2 $numSubdivJK2 -$y2 $z4 $y2 $z4 $y2 $z3 -$y2 $z3
	patch quad $IDconcCover $numSubdivIJ2 $numSubdivJK2 -$y2 -$z3 $y2 -$z3 $y2 -$z4 -$y2 -$z4
	patch quad $IDconcCover $numSubdivIJ3 $numSubdivJK3 $y3 -$z2 $y1 -$z2 $y1 $z2 $y3 $z2
	patch quad $IDconcCover $numSubdivIJ3 $numSubdivJK3 -$y1 -$z2 -$y3 -$z2 -$y3 $z2 -$y1 $z2
	patch quad $IDconcCover $numSubdivIJ4 $numSubdivJK4 $y2 -$z3 $y4 -$z3 $y4 $z3 $y2 $z3
	patch quad $IDconcCover $numSubdivIJ4 $numSubdivJK4 -$y4 -$z3 -$y2 -$z3 -$y2 $z3 -$y4 $z3
	#confined core
	patch quad $IDconcCore $numSubdivIJ5 $numSubdivJK5 -$y3 $z3 $y3 $z3 $y3 $z2 -$y3 $z2
	patch quad $IDconcCore $numSubdivIJ5 $numSubdivJK5 -$y3 -$z2 $y3 -$z2 $y3 -$z3 -$y3 -$z3
	patch quad $IDconcCore $numSubdivIJ6 $numSubdivJK6 $y4 -$z3 $y3 -$z3 $y3 $z3 $y4 $z3
	patch quad $IDconcCore $numSubdivIJ6 $numSubdivJK6 -$y3 -$z3 -$y4 -$z3 -$y4 $z3 -$y3 $z3
	#reinforcement
	#                                           $yStart $zStart $yEnd $zEnd
	layer straight $IDreinf $numFiber1 $areaFiber -$y3 $z2 $y3 $z2
	layer straight $IDreinf $numFiber1 $areaFiber -$y3 -$z2 $y3 -$z2
	layer straight $IDreinf $numFiber1 $areaFiber -$y3 $z3 $y3 $z3
	layer straight $IDreinf $numFiber1 $areaFiber -$y3 -$z3 $y3 -$z3
	layer straight $IDreinf $numFiber2 $areaFiber $y3 -$z4s $y3 $z4s
	layer straight $IDreinf $numFiber2 $areaFiber -$y3 -$z4s -$y3 $z4s
	layer straight $IDreinf $numFiber2 $areaFiber $y4 -$z4s $y4 $z4s
	layer straight $IDreinf $numFiber2 $areaFiber -$y4 -$z4s -$y4 $z4s
}

# assign torsional Stiffness for 3D Model
set SecTagTorsion 99;		# ID tag for torsional section behavior
set SecTag3D 3;			# ID tag for combined behavior for 3D model
uniaxialMaterial Elastic $SecTagTorsion $Ubig;	# define elastic torsional stiffness
section Aggregator $SecTag3D $SecTagTorsion T -section $SecTag;	# combine section properties

