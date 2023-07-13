## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: January 20, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# OCTAGONAL COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

################### Octagonal Column Fiber Cross Section Assignment ###################

# ColSecTag		      	-> Column's fiber element tag
# Dcol		      		-> Width of octagonal column (to flat sides)
# nLbar		      		-> Number of longitudinal bars
# DLbar					-> Diameter of longitudinal bars
# sTbar 				-> Spacing of transverse spiral reinforcement

proc BuildOctColSection {ColSecTag  Dcol  nLbar  DLbar  sTbar} {

	##### PRINT LINE TO CHECK #####
	# puts "Building octagonal column section for column $ColSecTag"
	
	if {$ColSecTag == 3010 || $ColSecTag == 5010 || $ColSecTag == 6010 || $ColSecTag == 8020 || $ColSecTag == 10010 || $ColSecTag == 11010 || $ColSecTag == 11020 || $ColSecTag == 12010 || $ColSecTag == 12020 || $ColSecTag == 12030 || $ColSecTag == 13010 || $ColSecTag == 14030} {
		set ftfact 1.0 
		##### PRINT LINE TO CHECK #####
		# puts "Setting tensile strength factor to $ftfact for column $ColSecTag"
		set fcfact 1.0
		##### PRINT LINE TO CHECK #####
		# puts "Setting compressive strength factor to $fcfact for column $ColSecTag"
	} else {
		set ftfact 1.0
		##### PRINT LINE TO CHECK #####
		# puts "Setting tensile strength factor to $ftfact for column $ColSecTag"
		set fcfact 1.0
		##### PRINT LINE TO CHECK #####
		# puts "Setting compressive strength factor to $fcfact for column $ColSecTag"
	};
	
	global kips in ksi_psi pi fce ec0 esp Ec Gc Es Esh fy fu esh esu; # Read in unit and constant variables as well as material property variables
	
	set IDconcCore  [expr $ColSecTag*10+1]
	set IDconcCover [expr $ColSecTag*10+2]
	set IDSteel     [expr $ColSecTag*10+3]
	set IDShear     [expr $ColSecTag*10+4]
	set IDTorsion   [expr $ColSecTag*10+5]
	
	# Column component dimensions
	set tcover 		[expr 2.0*$in];									# 2 inch cover width
	set Rcol		[expr $Dcol/2.0];								# Radius of octagonal column (to flat sides)
	set RcolDiag	[expr $Rcol/(cos($pi/8.0))];					# Radius of octagonal column (to corners)
	set Dcore		[expr $Dcol-2.0*$tcover];						# Diameter of circular core
	set Rcore		[expr $Dcore/2.0];								# Radius of circular core
	set Along		[expr $pi*$DLbar**2/4.0];						# Area of longitudinal reinforcement bar
	set DTbar		[expr 0.625*$in];								# Diameter of transverse spiral reinforcement bar (#5 Rebar)
	set Asp			[expr $pi*$DTbar**2/4.0];						# Area of transverse spiral reinforcement bar
	set Dtran		[expr $Dcore-$DTbar];							# Diameter of spiral of transverse spiral reinforcement
	set rho			[expr 4.0*$Asp/($Dtran*$sTbar)];				# Density of transverse spiral reinforcement
	set Dlong		[expr $Dcore-2*$DTbar-$DLbar];					# Diameter of ring of longitudinal reinforcement
	set Rlong		[expr $Dlong/2.0];								# Radius of ring of longitudinal reinforcement
	
	# Section Area Properties
	set Acol		[expr (2.0*$Dcol**2)/(1.0+sqrt(2.0))];			# Area of octagonal column section
	set Jcol		[expr 1.2762*$RcolDiag**4];	    				# Polar moment of inertia for octagonal column section
	set I3col		[expr 0.6381*$RcolDiag**4];	    				# Second moment of inertia, 1st transverse direction, for octagonal column section
	set I2col		[expr 0.6381*$RcolDiag**4];    					# Second moment of inertia, 2nd transverse direction, for octagonal column section
	
	# Material Definition - Steel (Steel02)
	set R0 18;									# control the transition from elastic to plastic branches
	set cR1 0.925;								# control the transition from elastic to plastic branches
	set cR2 0.15;								# control the transition from elastic to plastic branches
	uniaxialMaterial Steel02 $IDSteel $fy $Es 0.02 $R0 $cR1 $cR2
		
	# Material Definition - Concrete (Concrete02)
	# Compressive Properties	
	set ke 			[expr 1-$sTbar/$Dtran];					# Effective confinement strength coefficient from transverse reinforcement
	set f3e			[expr $ke*$rho*$fy/2]; 					# Effective confinement strength
	set fcc         [expr $fce*(-1.254+2.254*sqrt(1+7.94*$f3e/$fce)-2*$f3e/$fce)];	# Core compressive strength -compression
	set ecc         [expr $ec0*(1.0+5.0*($fcc/$fce-1.0))];	# Core strain at maximum strength -compression
	set ecu         [expr 0.004+$f3e/(4.0*$fce)];			# Core crushing (ultimate) strain -compression
	set lambda 0.1;											# Ratio between unloading slope at $eps2 and initial slope $Ec
	set xu 			[expr $ecu/$ecc];						# Mander equation parameter x
	set ru 			[expr $Ec/($Ec-$fcc/$ecc)]; 			# Mander equation parameter r
	set fcu         [expr $fcfact*$fcc*$xu*$ru/($ru-1+$xu**$ru)];	# Core crushing (ultimate) strength -compression
	set fcuu 		[expr 0.15*$fcc]; 						# Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
	set ecuu 		[expr $ecu+($ecu-$ecc)*($fcuu-$fcu)/($fcu-$fcc)]; # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression
	
	##### PRINT LINE TO CHECK #####
	# puts "f3e=$f3e, fcc=$fcc, ecc=$ecc, ecu=$ecu, fcu=$fcu, ecuu=$ecuu, fcuu=$fcuu";
	# puts "xu=$xu, ru=$ru, fcu=$fcu";
	
	# Tensile Properties -- zero tensile stiffness
	set ftU 0;	# Cover tensile strength +tension
	set ftC 0;	# Core tensile strength +tension
	set Ets 0;	# Tension softening stiffness
	# # Tensile Properties -- ACI tensile stiffness
	# set ftU [expr 7.5*sqrt($fce*$ksi_psi)/$ksi_psi];			# Cover tensile strength +tension
	# set ftC [expr $ftfact*7.5*sqrt($fcc*$ksi_psi)/$ksi_psi];	# Core tensile strength +tension
	# set Ets [expr $Ec/5.0];										# Tension softening stiffness
	
	##### PRINT LINE TO CHECK #####
	# puts "ftU=$ftU, ftC=$ftC, Ets=$Ets"
	
	# UniaxialMaterial Definition
	uniaxialMaterial Concrete02 $IDconcCover [expr -$fce] [expr -$ec0] [expr -0.1*$fce] [expr -$esp] $lambda $ftU $Ets;	# Cover concrete (unconfined)
	# puts "uniaxialMaterial Concrete02 $IDconcCover [expr -$fce] [expr -$ec0] [expr -0.1*$fce] [expr -$esp] $lambda $ftU $Ets;	# Cover concrete (unconfined)"
	
	##### PRINT LINE TO CHECK #####
	# puts "Cover fpc = [expr -$fce], epsc0 = [expr -$Efact*$ec0]"
	
	uniaxialMaterial Concrete02 $IDconcCore  [expr -$fcc] [expr -$ecc] [expr -$fcuu] [expr -$ecuu] $lambda $ftC $Ets;	# Core concrete (confined)
	# puts "uniaxialMaterial Concrete02 $IDconcCore  [expr -$fcc] [expr -$ecc] [expr -$fcuu] [expr -$ecuu] $lambda $ftC $Ets;	# Core concrete (confined)"
	
	##### PRINT LINE TO CHECK #####
	# puts "Core fpc = [expr -$fcc], epsc0 = [expr -$Efact*$ecc]"
	
	# Build Octagonal RC Column Section
	set numSubdivCirc 64; # Determines # of fibers in the circumferential direction, around the entire circumference
	set ColMatTag [expr 10*$ColSecTag];    # Column material cross-section tag.  Set equal to 10*Column's fiber element tag.
	uniaxialMaterial Elastic $IDTorsion [expr  0.2*$Gc*$Jcol];      # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1, fundamental period less than 0.7 sec, reduction is required *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED (this method of defining torsional stiffness is according to https://portwooddigital.com/2019/10/06/torsion-with-fiber-sections/).
	section Fiber $ColMatTag -torsion $IDTorsion {
		patch circ $IDconcCore $numSubdivCirc 5 0. 0. [expr 0.000e+00*$in] [expr $Rcore/2] 0.000e+00 3.600e+02; # Inner Core Patch, 5 radial fibers
		patch circ $IDconcCore $numSubdivCirc 10 0. 0. [expr $Rcore/2] [expr $Rcore] 0.000e+00 3.600e+02; # Outer Core Patch, 10 radial fibers
		set numSlices 8;	# Determines # of slices in each of the 8 sections of the octagon
		set numSubdivIJ 1;	# Determines # of fibers in the circumferential direction of the cover patch
		set numSubdivJK 4;	# Determines # of fibers in the radial direction of the cover patch
		for {set i 0} {$i < 8} {incr i 1} { # For each of the 8 sections of the octagon
			set startAngle	[expr $i*$pi/4+$pi/8]
			for {set j 0} {$j < $numSlices} {incr j 1} {
				set phi		[expr $pi/4/$numSlices]
				set sita1	[expr $startAngle + $j*$phi];	# Slice start angle
				set sita2	[expr $sita1 + $phi];			# Slice end angle
				set yI		[expr $Rcore*cos($sita2)]
				set zI		[expr $Rcore*sin($sita2)]
				set yJ		[expr $Rcore*cos($sita1)]
				set zJ		[expr $Rcore*sin($sita1)]
				set oR1		[expr $Rcol/cos($pi/8 - $j*$phi)]
				set oR2		[expr $Rcol/cos($pi/8 - ($j+1)*$phi)]
				set yK		[expr $oR1*cos($sita1)]
				set zK		[expr $oR1*sin($sita1)]
				set yL		[expr $oR2*cos($sita2)]
				set zL		[expr $oR2*sin($sita2)]
				patch quad $IDconcCover $numSubdivIJ $numSubdivJK $yI $zI $yJ $zJ $yK $zK $yL $zL; # Cover Patch connects the circular core to the octagonal cover
			}
		}
		layer circ $IDSteel $nLbar $Along 0. 0. $Rlong; # Longitudinal Bars
	}
	uniaxialMaterial Elastic $IDShear   [expr (9./10.)*$Gc*$Acol];	# Define elastic shear stiffness
	section Aggregator $ColSecTag $IDShear Vy $IDShear Vz -section $ColMatTag;
};