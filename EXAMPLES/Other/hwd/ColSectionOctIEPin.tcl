## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date:  January 20, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# OCTAGONAL COLUMN PIN SECTION DEFINITION PROCEDURE - *UNREINFORCED* INELASTIC PROPERTIES

################### Octagonal Column Fiber Cross Section Assignment ###################

# ColSecTag		      	-> Column's fiber element tag
# Dcol		      		-> Width of octagonal column (to flat sides)

proc BuildOctColPINSection {ColSecTag  Dcol} {

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
	
	set IDconc 		[expr $ColSecTag*10+2]
	set IDSteel     [expr $ColSecTag*10+3]
	set IDShear     [expr $ColSecTag*10+4]
	set IDTorsion   [expr $ColSecTag*10+5]
	
	# Column component dimensions
	set tcover 		[expr 2.0*$in];									# 2 inch cover width
	set Rcol		[expr $Dcol/2.0];								# Radius of octagonal column (to flat sides)
	set RcolDiag	[expr $Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0))];  	# Radius of octagonal column (to corners)
	set Dcore		[expr $Dcol-2.0*$tcover];						# Diameter of circular core
	set Rcore		[expr $Dcore/2.0];								# Radius of circular core

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

	# Concrete Properties
	set lambda 0.1;				# Ratio between unloading slope at $eps2 and initial slope $Ec
	set ftU 0;					# Tensile strength +tension
	set Ets 0;					# Tension softening stiffness
	
	##### PRINT LINE TO CHECK #####
	# puts "ftU=$ftU, Ets=$Ets"
	
	# UniaxialMaterial Definition
	uniaxialMaterial Concrete02 $IDconc [expr -$fce] [expr -$ec0] [expr -0.1*$fce] [expr -$esp] $lambda $ftU $Ets;	# Concrete (unconfined)
	# puts "uniaxialMaterial Concrete02 $IDconc [expr -$fce] [expr -$ec0] [expr -0.1*$fce] [expr -$esp] $lambda $ftU $Ets;	# Concrete (unconfined)"
	
	##### PRINT LINE TO CHECK #####
	# puts "fpc = [expr -$fce], epsc0 = [expr -$Efact*$ec0]"

	# Build Octagonal URC Column Section
	set numSubdivCirc 64; # Determines # of fibers in the circumferential direction, around the entire circumference
	set ColMatTag [expr 10*$ColSecTag];    # Column material cross-section tag.  Set equal to 10*Column's fiber element tag.
	uniaxialMaterial Elastic $IDTorsion [expr  0.2*$Gc*$Jcol];      # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1, fundamental period less than 0.7 sec, reduction is required *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED (this method of defining torsional stiffness is according to https://portwooddigital.com/2019/10/06/torsion-with-fiber-sections/).
	section Fiber $ColMatTag -torsion $IDTorsion {
		patch circ $IDconc $numSubdivCirc 5 0. 0. [expr 0.000e+00*$in] [expr $Rcore/2] 0.000e+00 3.600e+02; # Inner "Core" Patch, 5 radial fibers
		patch circ $IDconc $numSubdivCirc 10 0. 0. [expr $Rcore/2] [expr $Rcore] 0.000e+00 3.600e+02; # Outer "Core" Patch, 10 radial fibers
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
				patch quad $IDconc $numSubdivIJ $numSubdivJK $yI $zI $yJ $zJ $yK $zK $yL $zL; # Cover Patch connects the circular core to the octagonal cover
			}
		}
	layer circ $IDSteel 8 3.0 0. 0. 6.0; # Longitudinal Bars
	}
	uniaxialMaterial Elastic $IDShear   [expr (9./10.)*$Gc*$Acol];	# Define elastic shear stiffness
	section Aggregator $ColSecTag $IDShear Vy $IDShear Vz -section $ColMatTag;
};