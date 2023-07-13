# ----------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley, cchern@berkeley.edu.  Date: April 20, 2023.
# ----------------------------------------------------------------------------------
# MAIN MODEL CONFIGURATION SCRIPT
# Parameterized structural analysis model of the Hayward Hwy 580/238 Interchange bridge.
# To run with option files, type
	# "source hwd.tcl <path/to/conf_file.json> <path/to/analysis_file.json>
	# <path/to/results_requests_file> <path/to/output_directory> debug" into opensees
	# interpreter
# OR,
	# manual override with set override 1 and select options in the
	# "#+# HARD-CODED CONFIGURATION" section, then type "source hwd.tcl" into
	# opensees interpreter.
# This model is in units of kips, inches, and seconds.
# There are 3 dimensions and 6 degrees of freedom.
# 1/2 (X/Y) is the horizontal plane; 3 (Z) is the vertical axis
# global direction 1 (X): approximately longitudinal, east
# global direction 2 (Y): approximately transverse, north
# global direction 3 (Z): vertical, upwards
#
set override 1;  	# If you want to set the options manually, set this to 1 and
					# adjust the options in the "#+# HARD-CODED CONFIGURATION" section.
					# If this is set to 2, then this file is being sourced

# #+# SOME PROCEDURES --------------------------------------------------------------
proc py {args} {
    eval "[exec python {*}$args]"
}
proc write_modes {mode_file nmodes} {
  set fid_modes [open $mode_file w+]
  for {set m 1} {$m <= $nmodes} {incr m} {
    puts $fid_modes "$m:"
    foreach n [getNodeTags] {
      puts $fid_modes "  $n: \[[join [nodeEigenvector $n $m] {, }]\]";
    }
  }
  close $fid_modes
}
proc write_displacements {file_name {resp Disp}} {
  set fid [open "$file_name" "w+"]
  puts $fid "[getTime]:"
  foreach n [getNodeTags] {
	puts $fid "    $n: \[[join [node${resp} $n] {, }]\]";
  }
  close $fid;
}
proc make_column {tag inode jnode} {
    global column_linearity np Acol Ecol Gc Jcol I2col I3col
    if {$column_linearity == "nonlinear"} {
		set ele_args "$np $tag"
        set ele_cmd "nonlinearBeamColumn"
    } else {
		set ele_args "$Acol $Ecol $Gc $Jcol $I2col $I3col"
        set ele_cmd "elasticBeamColumn"
    }
    element $ele_cmd $tag $inode $jnode {*}$ele_args [expr int(floor($tag/100))]
}

if {$override == 0} {
#
# #+# PARSED CONFIGURATION ----------------------------------------------------------------
	set conf_file [lindex $argv 0]; 			# Set model configuration file
	set analysis_file [lindex $argv 1];			# Set analysis settings file
	set results_requests_file [lindex $argv 2]; # Set results requests file
	set output_directory [lindex $argv 3]; 		# Set data output directory
	file mkdir $output_directory;  				# Create data output directory if it doesn't exist
	if {[llength $argv] > 4} {
		set debug "--debug"; 					# Specify whether to print the specified options
	} else {
		set debug "";
	}
	eval [exec python -m CE58658.set_params $debug $conf_file]
	eval [exec python -m CE58658.set_params $debug $analysis_file]
	set damping_modes [split $damping_modes {,}];
	set damping_ratios [split $damping_ratios {,}];
	eval [exec python -m CE58658.set_params $debug $results_requests_file]

} elseif {$override == 2} {
	;

} else {
# #+# HARD-CODED CONFIGURATION ----------------------------------------------------------------
	set debug "--debug"; 												# print the specified options
	# model configuration
	set transformation Linear; 											# geometric transformation linearity for element deformations
	set column_linearity elastic; 										# column element material linearity
	set column_pins 3; 													# 1 - all rigid; 2 - all pinned; 3 - mixed rigid/pinned; 4 - all zerolength fiber sections; 5 - all integration point fiber sections
	set column_capbeam_joint none; 										# rigid offsets between tops of columns and column-cap beam joints
	set abutment_model linear; 											# abutment model ("none", "linear", "simplified", or "complex")
	set hinge_model linear; 											# in-span hinge model ("none", "linear", "simplified", or "complex")
	set Ec 3530.5; 														# concrete modulus of elasticity (ksi)
	set Ecol 3530.5; 													# column concrete modulus of elasticity (ksi)
	set Es 29000.0; 													# steel tensile modulus of elasticity (initial elastic tangent, ksi)
	set CGa 50.0; 														# abutment shear stiffness coefficient
	set CGh 60.0; 														# in-span hinge shear stiffness coefficient
	# analysis settings
	set dynamic_on 1; 													# turn the dynamic analysis on (1) or off (0)
	set dynamic_truncated 0; 											# truncate the dynamic analysis to the first t timesteps (1)
	set dynamic_timesteps 500;											# the t timesteps to which the the dynamic analysis is truncated
	set input_location "401"; 											# locations (node numbers) corresponding to ground motion input. 0 for multiple support excitation.
	set dynamic_integrator Newmark; 									# numerical integration method
	set damping_type rayleigh;        									# damping strategy
	set damping_modes "1,2"; 											# damping modes
	set damping_modes [split $damping_modes {,}];
	set damping_ratios "0.015,0.015"; 									# damping ratios
	set damping_ratios [split $damping_ratios {,}];
	set rayleigh_zerolength_on 1;      									# turn rayleigh damping on (1) for zerolength elements (abutment, hinge, and column pin springs)
	set dynamic_scale_factor 1.0;										# scale factor applied to input ground motion
	set record_zip "Records/58658_003_20210628_18.29.26.P_SanLo.zip"; 	# path to zip file containing the recorded motions for ground motion input
	set output_directory results;		 								# Set data output directory
	file mkdir $output_directory;  										# Create data output directory if it doesn't exist
	# results output
	set model_name "hwd_model"
	set modeling_matrix 1
	set runtime 1
	set eigen_modal_tracking 0
	set ssid_modal_tracking 0
	set compare_response_history 1
}
#
# #+# TIMER 0 ----------------------------------------------------------------------
set timers 	[open $output_directory/timers.json w];
set t0 [clock clicks -millisec];  	# Start timer for timing model and analyses definition and execution
#
# #+# INITIALIZE MODEL -------------------------------------------------------------
wipe;								# Clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6;	# Define the model builder, ndm=#dimension, ndf=#dofs
#
# #+# UNITS AND CONSTANTS ----------------------------------------------------------
# US UNITS
set kips  1.0;				# kips
set in    1.0;				# inches
set sec   1.0;				# seconds
# DEPENDENT UNITS
set ft      [expr $in*12];
set lb      [expr $kips/1000];
set ksi     [expr $kips/$in**2];  
set psi     [expr $lb/$in**2];
set ksi_psi [expr $ksi/$psi];
# CONSTANTS
set g		[expr 386.1*$in/$sec**2];
set pi      [expr acos(-1.0)];
#
# #+# MATERIAL PROPERTIES ----------------------------------------------------------
# CONCRETE PROPERTIES
set fc            [expr 3.5*$ksi];	                                # Default (general, if unspecified) 28-day concrete cylinder compressive strength   (+Tension, -Compression)
set fce           [expr 5.0*$ksi];									# Default (general, if unspecified) Unconfined compressive strength (max[5ksi or 1.3f'c])   (+Tension, -Compression)
set ec0           0.002;											# Unconfined strain at maximum strength (+Tension, -Compression)
set esp           0.005;											# Unconfined crushing (cover spalling) strain (+Tension, -Compression)
# set Ec            [expr 57000.0*sqrt($fce*$ksi_psi)/$ksi_psi];		# Default (general, if unspecified) Concrete Modulus of Elasticity
# set Ecol          [expr 57000.0*sqrt($fce*$ksi_psi)/$ksi_psi];		# Default (general, if unspecified) Concrete Modulus of Elasticity
# set Ecol          [expr 4395*$ksi];                         		# Average Column Concrete Core Test Modulus of Elasticity
set Ec 			  [expr $Ec*$ksi]
set Ecol 		  [expr $Ecol*$ksi]
set Uc            0.2;                                              # Poisson's ratio
set Gc            [expr $Ec/(2.0*(1.0+$Uc))];                       # Shear Modulus of Elasticity
set wconc	      [expr 143.96*$lb/$ft**3];                         # Normal Concrete Weight per Volume                              
set mconc	      [expr (143.96*$lb/$ft**3)/$g];                    # Normal Concrete Mass Weight per Volume
# REINFORCING STEEL PROPERTIES
# set Es           [expr 29000.0*$ksi];				# Steel Tensile Modulus of Elasticity (initial elastic tangent)
set Esh			 [expr 0.02*$Es];					# Tangent at initial strain hardening
set fy           [expr 68.0*$ksi];                	# Yield Strength
set fu			 [expr 95.0*$ksi];					# Ultimate strength
set esh          0.0075;                   			# Strain corresponding to initial strain hardening 
set esu          0.090;                    			# Strain at peak stress
#
# #+# NODES ------------------------------------------------------------------------
# See Haywardbridge_coordinates_NodeAsgn.xlsx for node coordinates and node assignment.
# COLUMN BENT NODES (Node numbering scheme: 100*Bent#+Node#)
node 201 11758.709 31062.549 [expr 140.750*$ft];
node 203 11758.709 31062.549 [expr 196.770*$ft];
node 204 11675.713 30944.873 [expr 197.490*$ft];
node 205 11571.969 30797.777 [expr 198.390*$ft];
node 207 11571.969 30797.777 [expr 140.750*$ft];
node 301 13326.879 30097.235 [expr 140.000*$ft];
node 303 13326.879 30097.235 [expr 199.570*$ft];
node 304 13227.933 29932.695 [expr 200.530*$ft];
node 305 12935.732 29446.785 [expr 203.370*$ft];
node 307 12935.732 29446.785 [expr 126.000*$ft];
node 401 15603.570 29883.836 [expr 137.000*$ft];
node 403 15603.570 29883.836 [expr 197.700*$ft];
node 404 15120.738 28901.035 [expr 203.145*$ft];
node 405 14983.164 28621.004 [expr 204.705*$ft];
node 407 14983.164 28621.004 [expr 119.000*$ft];
node 501 17372.788 28118.938 [expr 131.000*$ft];
node 503 17372.788 28118.938 [expr 204.004*$ft];
node 504 17308.901 27950.657 [expr 204.904*$ft];
node 505 17040.574 27243.878 [expr 208.684*$ft];
node 507 17040.574 27243.878 [expr 138.000*$ft];
node 601 19284.324 27591.751 [expr 142.000*$ft];
node 603 19284.324 27591.751 [expr 204.020*$ft];
node 604 19203.553 27315.309 [expr 205.460*$ft];
node 605 19047.059 26779.703 [expr 208.250*$ft];
node 607 19047.059 26779.703 [expr 144.000*$ft];
node 701 20759.780 27596.583 [expr 144.000*$ft];
node 703 20759.780 27596.583 [expr 202.334*$ft];
node 704 20609.812 26947.687 [expr 205.664*$ft];
node 705 20571.982 26784.002 [expr 206.504*$ft];
node 707 20571.982 26784.002 [expr 149.500*$ft];
node 801 22434.979 27020.686 [expr 145.000*$ft];
node 803 22434.979 27020.686 [expr 203.211*$ft];
node 804 22369.455 26605.829 [expr 205.311*$ft];
node 805 22332.013 26368.768 [expr 206.511*$ft];
node 807 22332.013 26368.768 [expr 154.000*$ft];
node 901 24104.380 26573.310 [expr 136.000*$ft];
node 903 24104.380 26573.310 [expr 203.028*$ft];
node 904 24088.125 26390.192 [expr 203.947*$ft];
node 905 24038.586 25832.224 [expr 206.748*$ft];
node 907 24038.586 25832.224 [expr 156.000*$ft];
node 1001 25885.617 26694.425 [expr 131.000*$ft];
node 1003 25885.617 26694.425 [expr 199.852*$ft];
node 1004 25877.420 26241.067 [expr 202.119*$ft];
node 1005 25872.332 25854.531 [expr 204.052*$ft];
node 1007 25872.332 25854.531 [expr 130.000*$ft];
node 1101 27396.787 26784.776 [expr 141.500*$ft];
node 1103 27396.787 26784.776 [expr 197.154*$ft];
node 1104 27421.866 26193.697 [expr 200.112*$ft];
node 1105 27432.903 25933.542 [expr 201.414*$ft];
node 1107 27432.903 25933.542 [expr 139.000*$ft];
node 1201 29283.281 26610.345 [expr 166.000*$ft];
node 1203 29283.281 26610.345 [expr 195.409*$ft];
node 1204 29300.043 26467.348 [expr 195.903*$ft];
node 1205 29326.939 26237.895 [expr 196.695*$ft];
node 1206 29372.201 25851.750 [expr 198.028*$ft];
node 1207 29370.596 25865.445 [expr 197.981*$ft];
node 1209 29370.596 25865.445 [expr 156.000*$ft];
node 1211 29326.939 26237.895 [expr 159.000*$ft];
node 1301 30813.126 26806.601 [expr 163.750*$ft];
node 1303 30813.126 26806.601 [expr 192.581*$ft];
node 1304 30749.330 26638.286 [expr 193.022*$ft];
node 1305 30685.534 26469.971 [expr 193.463*$ft];
node 1307 30685.534 26469.971 [expr 161.750*$ft];
node 1313 30439.057 25819.668 [expr 151.000*$ft];
node 1315 30439.057 25819.668 [expr 196.584*$ft];
node 1401 32390.923 27030.687 [expr 158.500*$ft];
node 1403 32390.923 27030.687 [expr 190.053*$ft];
node 1404 32222.481 26812.047 [expr 190.352*$ft];
node 1405 32054.040 26593.407 [expr 190.651*$ft];
node 1407 32054.040 26593.407 [expr 158.500*$ft];
node 1409 32222.481 26812.047 [expr 158.500*$ft];
node 1413 31434.869 25789.667 [expr 150.500*$ft];
node 1415 31434.869 25789.667 [expr 196.436*$ft];
# COLUMN JOINT OFFSET NODES
if {$column_capbeam_joint == "rigidlink"} {
	node 202 11758.709 31062.549 [expr 192.776*$ft];
	node 206 11571.969 30797.777 [expr 194.397*$ft];
	node 302 13326.879 30097.235 [expr 195.330*$ft];
	node 306 12935.732 29446.785 [expr 199.070*$ft];
	node 402 15603.570 29883.836 [expr 189.910*$ft];
	node 406 14983.164 28621.004 [expr 196.930*$ft];
	node 502 17372.788 28118.938 [expr 198.740*$ft];
	node 506 17040.574 27243.878 [expr 203.420*$ft];
	node 602 19284.324 27591.751 [expr 199.780*$ft];
	node 606 19047.059 26779.703 [expr 204.010*$ft];
	node 702 20759.780 27596.583 [expr 198.470*$ft];
	node 706 20571.982 26784.002 [expr 202.650*$ft];
	node 802 22434.979 27020.686 [expr 199.690*$ft];
	node 806 22332.013 26368.768 [expr 202.990*$ft];
	node 902 24104.380 26573.310 [expr 199.000*$ft];
	node 906 24038.586 25832.224 [expr 202.720*$ft];
	node 1002 25885.617 26694.425 [expr 194.810*$ft];
	node 1006 25872.332 25854.531 [expr 199.010*$ft];
	node 1102 27396.787 26784.776 [expr 192.090*$ft];
	node 1106 27432.903 25933.542 [expr 196.350*$ft];
	node 1202 29283.281 26610.345 [expr 192.159*$ft];
	node 1208 29370.596 25865.445 [expr 193.445*$ft];
	node 1210 29326.939 26237.895 [expr 194.731*$ft];
	node 1302 30813.126 26806.601 [expr 189.331*$ft];
	node 1306 30685.534 26469.971 [expr 190.215*$ft];
	node 1314 30439.057 25819.668 [expr 193.333*$ft];
	node 1402 32390.923 27030.687 [expr 186.803*$ft];
	node 1406 32054.040 26593.407 [expr 187.402*$ft];
	node 1408 32222.481 26812.047 [expr 187.102*$ft];
	node 1414 31434.869 25789.667 [expr 193.182*$ft];
}
# COLUMN PIN NODES
# (Node numbering scheme: 1000*Bent#+10*Node#,
# i.e. add 0 to the end of column bent node number)
if {$column_pins != 1 && $column_pins !=5} {
	if {$column_capbeam_joint == "rigidlink"} {
		node 2010 11758.709 31062.549 [expr 140.750*$ft];
		node 2070 11571.969 30797.777 [expr 140.750*$ft];
		node 3020 13326.879 30097.235 [expr 195.330*$ft];
		node 3060 12935.732 29446.785 [expr 199.070*$ft];
		node 4020 15603.570 29883.836 [expr 189.910*$ft];
		node 4060 14983.164 28621.004 [expr 196.930*$ft];
		node 5020 17372.788 28118.938 [expr 198.740*$ft];
		node 5060 17040.574 27243.878 [expr 203.420*$ft];
		node 6020 19284.324 27591.751 [expr 199.780*$ft];
		node 6060 19047.059 26779.703 [expr 204.010*$ft];
		node 7020 20759.780 27596.583 [expr 198.470*$ft];
		node 7060 20571.982 26784.002 [expr 202.650*$ft];
		node 8020 22434.979 27020.686 [expr 199.690*$ft];
		node 8060 22332.013 26368.768 [expr 202.990*$ft];
		node 9020 24104.380 26573.310 [expr 199.000*$ft];
		node 9060 24038.586 25832.224 [expr 202.720*$ft];
		node 10020 25885.617 26694.425 [expr 194.810*$ft];
		node 10060 25872.332 25854.531 [expr 199.010*$ft];
		node 11020 27396.787 26784.776 [expr 192.090*$ft];
		node 11060 27432.903 25933.542 [expr 196.350*$ft];
		node 12010 29283.281 26610.345 [expr 166.000*$ft];
		node 12090 29370.596 25865.445 [expr 156.000*$ft];
		node 12110 29326.939 26237.895 [expr 159.000*$ft];
		node 13010 30813.126 26806.601 [expr 163.750*$ft];
		node 13070 30685.534 26469.971 [expr 161.750*$ft];
		node 14010 32390.923 27030.687 [expr 158.500*$ft];
		node 14070 32054.040 26593.407 [expr 158.500*$ft];
		node 14090 32222.481 26812.047 [expr 158.500*$ft];
	} else {
		node 2010 11758.709 31062.549 [expr 140.750*$ft];
		node 2070 11571.969 30797.777 [expr 140.750*$ft];
		node 3030 13326.879 30097.235 [expr 199.570*$ft];
		node 3050 12935.732 29446.785 [expr 203.370*$ft];
		node 4030 15603.570 29883.836 [expr 197.700*$ft];
		node 4050 14983.164 28621.004 [expr 204.705*$ft];
		node 5030 17372.788 28118.938 [expr 204.004*$ft];
		node 5050 17040.574 27243.878 [expr 208.684*$ft];
		node 6030 19284.324 27591.751 [expr 204.020*$ft];
		node 6050 19047.059 26779.703 [expr 208.250*$ft];
		node 7030 20759.780 27596.583 [expr 202.334*$ft];
		node 7050 20571.982 26784.002 [expr 206.504*$ft];
		node 8030 22434.979 27020.686 [expr 203.211*$ft];
		node 8050 22332.013 26368.768 [expr 206.511*$ft];
		node 9030 24104.380 26573.310 [expr 203.028*$ft];
		node 9050 24038.586 25832.224 [expr 206.748*$ft];
		node 10030 25885.617 26694.425 [expr 199.852*$ft];
		node 10050 25872.332 25854.531 [expr 204.052*$ft];
		node 11030 27396.787 26784.776 [expr 197.154*$ft];
		node 11050 27432.903 25933.542 [expr 201.414*$ft];
		node 12010 29283.281 26610.345 [expr 166.000*$ft];
		node 12090 29370.596 25865.445 [expr 156.000*$ft];
		node 12110 29326.939 26237.895 [expr 159.000*$ft];
		node 13010 30813.126 26806.601 [expr 163.750*$ft];
		node 13070 30685.534 26469.971 [expr 161.750*$ft];
		node 14010 32390.923 27030.687 [expr 158.500*$ft];
		node 14070 32054.040 26593.407 [expr 158.500*$ft];
		node 14090 32222.481 26812.047 [expr 158.500*$ft];
	}
}
# DECK NODES
# (Node numbering scheme: 10000*1stBent#+Node#)
node 10001 10956.246 31471.424 [expr 195.779*$ft];
node 10002 11134.376 31337.427 [expr 196.224*$ft];
node 10003 11313.673 31204.996 [expr 196.658*$ft];
node 10004 11494.123 31074.141 [expr 197.082*$ft];
node 20001 11980.192 30733.432 [expr 198.159*$ft];
node 20002 12287.718 30526.451 [expr 198.794*$ft];
node 20003 12598.228 30323.972 [expr 199.401*$ft];
node 20004 12911.655 30126.039 [expr 199.978*$ft];
node 30001 13599.389 29713.586 [expr 201.129*$ft];
node 30002 13974.506 29500.805 [expr 201.691*$ft];
node 30003 14353.175 29294.414 [expr 202.215*$ft];
node 30004 14735.289 29094.471 [expr 202.700*$ft];
node 40001 15551.083 28694.603 [expr 203.593*$ft];
node 40002 15985.225 28496.280 [expr 203.992*$ft];
node 40003 16423.012 28306.138 [expr 204.344*$ft];
node 40004 16864.289 28124.242 [expr 204.648*$ft];
node 50001 17604.098 27840.823 [expr 205.081*$ft];
node 50002 18060.724 27678.668 [expr 205.225*$ft];
node 50003 18439.748 27551.567 [expr 205.335*$ft];
node 50004 18820.722 27430.437 [expr 205.412*$ft];
node 60001 19483.087 27235.365 [expr 205.533*$ft];
node 60002 19763.517 27158.622 [expr 205.592*$ft];
node 60003 20044.806 27085.088 [expr 205.634*$ft];
node 60004 20326.916 27014.773 [expr 205.658*$ft];
node 70001 20959.743 26869.407 [expr 205.648*$ft];
node 70002 21310.743 26796.068 [expr 205.604*$ft];
node 70003 21662.742 26727.682 [expr 205.534*$ft];
node 70004 22015.669 26664.265 [expr 205.436*$ft];
node 80001 22712.048 26554.108 [expr 205.085*$ft];
node 80002 23122.231 26498.529 [expr 204.834*$ft];
node 80003 23399.187 26464.685 [expr 204.557*$ft];
node 80004 23743.626 26427.172 [expr 204.258*$ft];
node 90001 24446.148 26360.362 [expr 203.584*$ft];
node 90002 24804.172 26330.532 [expr 203.220*$ft];
node 90003 25162.197 26300.702 [expr 202.854*$ft];
node 90004 25520.222 26270.872 [expr 202.486*$ft];
node 100001 26186.309 26231.593 [expr 201.719*$ft];
node 100002 26495.198 26222.119 [expr 201.319*$ft];
node 100003 26804.087 26212.645 [expr 200.917*$ft];
node 100004 27112.976 26203.171 [expr 200.515*$ft];
node 110001 27802.880 26202.537 [expr 199.488*$ft];
node 110002 28183.895 26211.376 [expr 198.835*$ft];
node 110003 28564.909 26220.216 [expr 198.143*$ft];
node 110004 28945.924 26229.055 [expr 197.397*$ft];
node 120001 29526.646 26494.061 [expr 195.445*$ft];
node 120002 29879.758 26535.723 [expr 194.631*$ft];
node 120003 30169.615 26569.911 [expr 193.995*$ft];
node 120004 30459.472 26604.098 [expr 193.360*$ft];
node 120005 29603.224 25845.334 [expr 197.570*$ft];
node 120006 29798.944 25838.917 [expr 197.342*$ft];
node 120007 30012.315 25832.501 [expr 197.051*$ft];
node 120008 30225.686 25826.084 [expr 196.808*$ft];
node 130001 31043.962 26673.037 [expr 192.337*$ft];
node 130002 31338.594 26707.787 [expr 191.745*$ft];
node 130003 31633.226 26742.538 [expr 191.188*$ft];
node 130004 31927.858 26777.289 [expr 190.665*$ft];
node 130005 30638.228 25813.678 [expr 196.460*$ft];
node 130006 30837.398 25807.689 [expr 196.384*$ft];
node 130007 31036.569 25801.699 [expr 196.353*$ft];
node 130008 31235.740 25795.709 [expr 196.368*$ft];
node 140001 32552.005 26850.906 [expr 189.683*$ft];
node 140002 32881.519 26889.773 [expr 189.240*$ft];
node 140003 33211.033 26928.640 [expr 188.868*$ft];
node 140004 33540.547 26967.506 [expr 188.566*$ft];
node 140005 31783.829 25779.163 [expr 196.639*$ft];
node 140006 32132.789 25768.658 [expr 196.984*$ft];
node 140007 32481.748 25758.154 [expr 197.462*$ft];
node 140008 32830.708 25747.650 [expr 198.074*$ft];
# IN-SPAN HINGE NODES
# (Node numbering scheme: 100000*1stBent#+10*Node#,
# i.e. add 0 to the end of deck node number)
if {$hinge_model != "none"} {
	node 500010 17604.098 27840.823 [expr 205.081*$ft];
	node 800040 23743.626 26427.172 [expr 204.258*$ft];
	node 1200010 29526.646 26494.061 [expr 195.445*$ft];
	node 1200050 29603.224 25845.334 [expr 197.570*$ft];
}
# ABUTMENT NODES (Node numbering scheme: 100*Bent#+Node#)
node 1010 10779.297 31606.976 [expr 194.965*$ft];
node 1020 10935.542 31809.096 [expr 193.626*$ft];
node 1030 10623.052 31404.856 [expr 196.183*$ft];
node 15010 33870.061 27006.373 [expr 188.346*$ft];
node 15020 34094.051 27297.117 [expr 188.604*$ft];
node 15030 33646.070 26715.629 [expr 188.089*$ft];
node 15040 33179.709 25737.199 [expr 198.822*$ft];
node 15050 33257.634 25838.348 [expr 198.650*$ft];
node 15060 33101.704 25635.948 [expr 198.994*$ft];
# ABUTMENT SPRING NODES, FIXED ENDS
# (Node numbering scheme: Node numbering scheme: 100*Bent#+Node#+1)
if {$abutment_model != "none"} {
	node 1021 10935.542 31809.096 [expr 193.626*$ft];
	node 1031 10623.052 31404.856 [expr 196.183*$ft];
	node 15021 34094.051 27297.117 [expr 188.604*$ft];
	node 15031 33646.070 26715.629 [expr 188.089*$ft];
	node 15051 33257.634 25838.348 [expr 198.650*$ft];
	node 15061 33101.704 25635.948 [expr 198.994*$ft];
}
# #+# CONSTRAINTS ------------------------------------------------------------------
# Fixed base for the columns
fix 201  1 1 1 1 1 1;
fix 207  1 1 1 1 1 1;
fix 301  1 1 1 1 1 1;
fix 307  1 1 1 1 1 1;
fix 401  1 1 1 1 1 1;
fix 407  1 1 1 1 1 1;
fix 501  1 1 1 1 1 1;
fix 507  1 1 1 1 1 1;
fix 601  1 1 1 1 1 1;
fix 607  1 1 1 1 1 1;
fix 701  1 1 1 1 1 1;
fix 707  1 1 1 1 1 1;
fix 801  1 1 1 1 1 1;
fix 807  1 1 1 1 1 1;
fix 901  1 1 1 1 1 1;
fix 907  1 1 1 1 1 1;
fix 1001  1 1 1 1 1 1;
fix 1007  1 1 1 1 1 1;
fix 1101  1 1 1 1 1 1;
fix 1107  1 1 1 1 1 1;
fix 1201  1 1 1 1 1 1;
fix 1209  1 1 1 1 1 1;
fix 1211  1 1 1 1 1 1;
fix 1301  1 1 1 1 1 1;
fix 1307  1 1 1 1 1 1;
fix 1313  1 1 1 1 1 1;
fix 1401  1 1 1 1 1 1;
fix 1407  1 1 1 1 1 1;
fix 1409  1 1 1 1 1 1;
fix 1413  1 1 1 1 1 1;
# Fixed abutment spring ends
if {$abutment_model != "none"} {
	fix 1021  1 1 1 1 1 1;
	fix 1031  1 1 1 1 1 1;
	fix 15021 1 1 1 1 1 1;
	fix 15031 1 1 1 1 1 1;
	fix 15051 1 1 1 1 1 1;
	fix 15061 1 1 1 1 1 1;
} else {
	fix 1020  1 1 1 1 1 1;
	fix 1030  1 1 1 1 1 1;
	fix 15020 1 1 1 1 1 1;
	fix 15030 1 1 1 1 1 1;
	fix 15050 1 1 1 1 1 1;
	fix 15060 1 1 1 1 1 1;
}
#
# #+# COLUMN FIBER CROSS-SECTION ASSIGNMENT ----------------------------------------
# (Column Numbering scheme: 1000*Bent#+10*Column#)
# (Column# (Facing East): 1=Left(North), 2=Right(South), 3=Center, 4=Single)
# Column fiber element tag (ColSecTag) follows column numbering scheme.
# Read in the column heights:
set fpHcol [open "./Dimensions/Hcol.txt" r];
set HcolList [read $fpHcol];
close $fpHcol
# Read in the number of longitudinal bars for each column:
set fpnLbar [open "./Dimensions/nLbar.txt" r];
set nLbarList [read $fpnLbar];
close $fpnLbar
# Read in the diameter of longitudinal bars for each column:
set fpDLbar [open "./Dimensions/DLbar.txt" r];
set DLbarList [read $fpDLbar];
close $fpDLbar
# Read in the spacing of transverse spiral reinforcement for each column:
set fpsTbar [open "./Dimensions/sTbar.txt" r];
set sTbarList [read $fpsTbar];
close $fpsTbar
# Procedure for column fiber cross-section assignment (INELASTIC material properties):
source ColSectionOctIEu.tcl
source ColSectionOctWideIEu.tcl
# BUILD COLUMN SECTIONS
if {$column_linearity == "nonlinear"} {
# Bents 2-11
for {set ib 2} {$ib <= 11} {incr ib 1} {
	set Dcol [expr 84.0*$in];
	set ic 1;											# North column
	set ColSecTag [expr 1000*$ib+10*$ic];				# Column's fiber element tag. Follows column numbering scheme.
	set Hcol [lindex $HcolList [expr ($ib-2)*2]];
	set nLbar [lindex $nLbarList [expr ($ib-2)*2]];
	set DLbar [lindex $DLbarList [expr ($ib-2)*2]];
	set sTbar [lindex $sTbarList [expr ($ib-2)*2]];
	BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
	set ic 2;											# South column
	set ColSecTag [expr 1000*$ib+10*$ic];				# Column's fiber element tag. Follows column numbering scheme.
	set Hcol [lindex $HcolList [expr ($ib-2)*2+1]];
	set nLbar [lindex $nLbarList [expr ($ib-2)*2+1]];
	set DLbar [lindex $DLbarList [expr ($ib-2)*2+1]];
	set sTbar [lindex $sTbarList [expr ($ib-2)*2+1]];
	BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
}
# Bent 12
set ib 12;
set Dcol [expr 66.0*$in];
set ic 1;												# North column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 20];
set nLbar [lindex $nLbarList 20];
set DLbar [lindex $DLbarList 20];
set sTbar [lindex $sTbarList 20];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
set ic 2;												# South column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 22];
set nLbar [lindex $nLbarList 22];
set DLbar [lindex $DLbarList 22];
set sTbar [lindex $sTbarList 22];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
set ic 3;												# Center Column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 21];
set nLbar [lindex $nLbarList 21];
set DLbar [lindex $DLbarList 21];
set sTbar [lindex $sTbarList 21];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
# Bent 13 NE
set ib 13;
set Dcol [expr 48.0*$in];
set ic 1;												# North column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 23];
set nLbar [lindex $nLbarList 23];
set DLbar [lindex $DLbarList 23];
set sTbar [lindex $sTbarList 23];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
set ic 2;												# South column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 24];
set nLbar [lindex $nLbarList 24];
set DLbar [lindex $DLbarList 24];
set sTbar [lindex $sTbarList 24];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
set ib 13;
set Dcol [expr 48.0*$in]; 								# Shorter Width of octagonal column (to flat sides)
set ic 4;												# Single Column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 25]; 						# Column Height
set nLbar [lindex $nLbarList 25]; 						# Number of main (outer) longitudinal bars
set DLbar [lindex $DLbarList 25]; 						# Diameter of main (outer) longitudinal bars
set sTbar [lindex $sTbarList 25]; 						# Spacing of transverse spiral reinforcement
set Wcol [expr 72.0*$in]; 								# Longer Width of octagonal column (to flat sides)
set nLbar2 8; 											# Number of secondary (inner) longitudinal bars
set DLbar2 [expr 0.625*$in]; 							# Diameter of secondary (inner) longitudinal bars (#5 rebar)
BuildWideOctColSection $ColSecTag  $Dcol  $Wcol  $nLbar  $nLbar2  $DLbar  $DLbar2  $sTbar; # Fiber cross-section
# Bent 14 NE
set ib 14;
set Dcol [expr 48.0*$in];
set ic 1;												# North column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 26];
set nLbar [lindex $nLbarList 26];
set DLbar [lindex $DLbarList 26];
set sTbar [lindex $sTbarList 26];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
set ic 2;												# South column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 28];
set nLbar [lindex $nLbarList 28];
set DLbar [lindex $DLbarList 28];
set sTbar [lindex $sTbarList 28];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
set ic 3;												# Center Column, facing east
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 27];
set nLbar [lindex $nLbarList 27];
set DLbar [lindex $DLbarList 27];
set sTbar [lindex $sTbarList 27];
BuildOctColSection $ColSecTag  $Dcol  $nLbar  $DLbar $sTbar; # Fiber cross-section
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
set ib 14;
set Dcol [expr 48.0*$in]; 								# Shorter Width of octagonal column (to flat sides)
set ic 4; 												# Single Column
set ColSecTag [expr 1000*$ib+10*$ic];					# Column's fiber element tag. Follows column numbering scheme.
set Hcol [lindex $HcolList 29]; 						# Column Height
set nLbar [lindex $nLbarList 29]; 						# Number of main (outer) longitudinal bars
set DLbar [lindex $DLbarList 29]; 						# Diameter of main (outer) longitudinal bars
set sTbar [lindex $sTbarList 29]; 						# Spacing of transverse spiral reinforcement
set Wcol [expr 72.0*$in]; 								# Longer Width of octagonal column (to flat sides)
set nLbar2 8; 											# Number of secondary (inner) longitudinal bars
set DLbar2 [expr 0.625*$in]; 							# Diameter of secondary (inner) longitudinal bars (#5 rebar)
BuildWideOctColSection $ColSecTag  $Dcol  $Wcol  $nLbar  $nLbar2  $DLbar  $DLbar2  $sTbar; # Fiber cross-section
}
#
# #+# COLUMN **PIN** FIBER CROSS-SECTION ASSIGNMENT ---------------------------------------------------------------------
#
if {$column_pins == 5} {
	source ColSectionOctIEPin.tcl
	set Dcol [expr 84.0*$in];
	BuildOctColPINSection 290000  $Dcol; # Fiber cross-section for columns at Bents 2-11 (each have diameter 84.0 inches)
	set Dcol [expr 66.0*$in];
	BuildOctColPINSection 1290000  $Dcol; # Fiber cross-section for columns at Bent 12 (each have diameter 66.0 inches) 
	set Dcol [expr 48.0*$in];
	BuildOctColPINSection 1390000  $Dcol; # Fiber cross-section for columns at Bents 13-14 (each have diameter 48.0 inches)
}
#
# #+# ELEMENTS ---------------------------------------------------------------------
# COLUMN ELEMENTS
# (Numbering scheme: Bottom Length, Hinge: 1000*Bent#+10*Column#;  Top Length, Rigid: 1000*Bent#+10*Column#+1)
# (Column# (Facing East): 1=Left(North), 2=Right(South), 3=Center, 4=Single)
# Column section tag (ColSecTag) follows column bent numbering scheme.
# Read in the local axes transformation vectors (vecxz) for the columns.
# vecxz is perpendicular to the length of cap beam, pointing approximately east.
# this makes the local y (direction 1) axis of the column section parallel to the cap beam, pointing approximately south, and
# this makes the local z (direction 2) axis of the column section perpendicular to the cap beam, pointing approximately east.
set fpvecxzX [open "./Dimensions/vecxzXcol.txt" r];
set vecxzXList [read $fpvecxzX];
close $fpvecxzX
set fpvecxzY [open "./Dimensions/vecxzYcol.txt" r];
set vecxzYList [read $fpvecxzY];
close $fpvecxzY
set np 4;														# number of Gauss integration points for nonlinear curvature distribution-- np=2 for linear distribution ok
set locations 0.0; 												# location of integrated pin fiber section
for {set ip 1} {$ip < $np} {incr ip 1} {
	lappend locations [expr $ip/(1.0*$np-1)];}; # locations of Fixed Location integration points for nonlinear curvature distribution
# Bent 2
set Dcol        [expr 84.0*$in];								# Width of octagonal column (to flat sides)
set RcolDiag	[expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];	# Radius of octagonal column (to corners)
set Acol		[expr (2.0*$Dcol**2)/(1.0+sqrt(2.0))];			# Area of octagonal column section
set Jcol		[expr 1.2762*$RcolDiag**4];	    				# Polar moment of inertia for octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    				# Second moment of inertia, 1st transverse direction, for octagonal column section
set I3col		[expr 0.6381*$RcolDiag**4];    					# Second moment of inertia, 2nd transverse direction, for octagonal column section
set vecxzX 		[lindex $vecxzXList [expr 0]]; 					# X component of vecxz vector for local axis transformation
set vecxzY 		[lindex $vecxzYList [expr 0]]; 					# Y component of vecxz vector for local axis transformation
geomTransf $transformation 20 $vecxzX $vecxzY 0; 				# vecxz is perpendicular to the length of the cap beam.
if {$column_capbeam_joint == "rigidlink"} {
	set topN 202
	set topS 206
	rigidLink beam 202 203; 
	rigidLink beam 206 205; 	
} else {
	set topN 203
	set topS 205
}
if {$column_pins != 1 && $column_pins !=5} {
    make_column 2010 2010 $topN; 
    make_column 2020 2070 $topS; 
} elseif {$column_pins == 1} {
    make_column 2010 201 $topN; 
    make_column 2020 207 $topS; 
} elseif {$column_pins == 5} {
	set secTags 290000
	for {set ip 1} {$ip < $np} {incr ip 1} {lappend secTags 2010;};
	set integration "FixedLocation $np $secTags $locations"
	element forceBeamColumn  	2010   201   $topN    20  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
	set secTags 290000
	for {set ip 1} {$ip < $np} {incr ip 1} {
		lappend secTags 2020;};
	set integration "FixedLocation $np $secTags $locations"
	element forceBeamColumn 	2020   207   $topS    20  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
}
# Bents 3-11
for {set ib 3} {$ib <= 11} {incr ib 1} {
	set vecxzX 		[lindex $vecxzXList [expr $ib-2]]; 					# X component of vecxz vector for local axis transformation
	set vecxzY 		[lindex $vecxzYList [expr $ib-2]]; 					# Y component of vecxz vector for local axis transformation
	geomTransf $transformation [expr 10*$ib] $vecxzX $vecxzY 0; 		# vecxz is perpendicular to the length of the cap beam.
	if {$column_capbeam_joint == "rigidlink"} {
		set topN 2
		set topS 6
		if {$column_pins != 1 && $column_pins !=5} {
			rigidLink beam [expr 1000*$ib+20] [expr 100*$ib+3];	
			rigidLink beam [expr 1000*$ib+60] [expr 100*$ib+5];	
		} else {
			rigidLink beam [expr 100*$ib+2] [expr 100*$ib+3]; 
			rigidLink beam [expr 100*$ib+6] [expr 100*$ib+5]; 
		}  
	} else {
		set topN 3
		set topS 5
	}
    if {$column_pins !=5} {
		make_column [expr 1000*$ib+10] [expr 100*$ib+1] [expr 100*$ib+$topN]; 
		make_column [expr 1000*$ib+20] [expr 100*$ib+7] [expr 100*$ib+$topS]; 
	} else {
		set secTags 290000
		for {set ip 1} {$ip < $np} {incr ip 1} {
			lappend secTags [expr 1000*$ib+10];};
		set integration "FixedLocation $np [lreverse $secTags] $locations"
		element forceBeamColumn  	[expr 1000*$ib+10] [expr 100*$ib+1] [expr 100*$ib+$topN] [expr 10*$ib]  $integration; 	# Fixed Location integration points, point at TOP has PIN fiber section)
		set secTags 290000
		for {set ip 1} {$ip < $np} {incr ip 1} {
			lappend secTags [expr 1000*$ib+20];};
		set integration "FixedLocation $np [lreverse $secTags] $locations"
		element forceBeamColumn  	[expr 1000*$ib+20] [expr 100*$ib+7] [expr 100*$ib+$topS] [expr 10*$ib]  $integration; 	# Fixed Location integration points, point at TOP has PIN fiber section)	
	} 
}
# Bent 12
set Dcol        [expr 66.0*$in];								# Width of octagonal column (to flat sides)
set RcolDiag	[expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];	# Radius of octagonal column (to corners)
set Acol		[expr (2.0*$Dcol**2)/(1.0+sqrt(2.0))];			# Area of octagonal column section
set Jcol		[expr 1.2762*$RcolDiag**4];	    				# Polar moment of inertia for octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    				# Second moment of inertia, 1st transverse direction, for octagonal column section
set I3col		[expr 0.6381*$RcolDiag**4];    					# Second moment of inertia, 2nd transverse direction, for octagonal column section
set vecxzX 		[lindex $vecxzXList 10]; 						# X component of vecxz vector for local axis transformation
set vecxzY 		[lindex $vecxzYList 10]; 						# Y component of vecxz vector for local axis transformation
geomTransf $transformation 120 $vecxzX $vecxzY 0; 				# vecxz is perpendicular to the length of the cap beam.
if {$column_capbeam_joint == "rigidlink"} {
	set topN 1202
	set topS 1208
	set topM 1210
	rigidLink beam 1202 1203; 
	rigidLink beam 1208 1207; 
	rigidLink beam 1210 1205; 	
} else {
	set topN 1203
	set topS 1207
	set topM 1205
}
if {$column_pins != 1 && $column_pins !=5} {
    make_column 12010 12010 $topN; 
    make_column 12020 12090 $topS; 
    make_column 12030 12110 $topM; 
} elseif {$column_pins == 1} {
    make_column 12010 1201 $topN; 
    make_column 12020 1209 $topS; 
    make_column 12030 1211 $topM; 
} elseif {$column_pins == 5} {
	set secTags 1290000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 12010;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	12010   1201   $topN    120  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
    set secTags 1290000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 12020;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	12020   1209   $topS    120  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
    set secTags 1290000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 12030;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	12030   1211   $topM    120  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
}
# Bent 13 NE
set Dcol        [expr 48.0*$in];								# Width of octagonal column (to flat sides)
set RcolDiag	[expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];	# Radius of octagonal column (to corners)
set Acol		[expr (2.0*$Dcol**2)/(1.0+sqrt(2.0))];			# Area of octagonal column section
set Jcol		[expr 1.2762*$RcolDiag**4];	    				# Polar moment of inertia for octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    				# Second moment of inertia, 1st transverse direction, for octagonal column section
set I3col		[expr 0.6381*$RcolDiag**4];    					# Second moment of inertia, 2nd transverse direction, for octagonal column section
set vecxzX 		[lindex $vecxzXList 11]; 						# X component of vecxz vector for local axis transformation
set vecxzY 		[lindex $vecxzYList 11]; 						# Y component of vecxz vector for local axis transformation
geomTransf $transformation 130 $vecxzX $vecxzY 0; 				# vecxz is perpendicular to the length of the cap beam.
if {$column_capbeam_joint == "rigidlink"} {
	set topN 1302
	set topS 1306
	rigidLink beam 1302 1303; 
	rigidLink beam 1306 1305; 	
} else {
	set topN 1303
	set topS 1305
}
if {$column_pins != 1 && $column_pins !=5} {
    make_column 13010 13010 $topN;
    make_column 13020 13070 $topS;
} elseif {$column_pins == 1} {
    make_column 13010 1301 $topN;
    make_column 13020 1307 $topS;
} elseif {$column_pins == 5} {
    set secTags 1390000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 13010;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	13010   1301   $topN    130  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
    set secTags 1390000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 13020;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	13020   1307   $topS    130  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
}
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
set Wcol        [expr 72.0*$in]; 										# Longer Width of octagonal column (to flat sides)
set spO         [expr ($Wcol-$Dcol)/2.0]; 								# Offset of octagonal sections from centroid (along horizontal axis)
set Acol		[expr (2*$Dcol**2)/(1+sqrt(2)) + ($Wcol-$Dcol)*$Dcol]; 	# Area of WIDE octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    						# Second moment of inertia about the 1st transverse (local y, horizontal) axis, for WIDE octagonal column section
set I3col		[expr $I2col+$Acol*$spO**2+($Wcol-$Dcol)*$Dcol**3/12]; 	# Second moment of inertia about the 2nd transverse (local z, vertical) axis, for WIDE octagonal column section
set Jcol		[expr $I2col+$I3col];		    						# Polar moment of inertia for WIDE octagonal column section
if {$column_capbeam_joint == "rigidlink"} {
	make_column 13040 1313 1314;
	rigidLink beam 1314 1315;	
} else {
	make_column 13040 1313 1315;
}
# Bent 14 NE
set Dcol        [expr 48.0*$in];								# Width of octagonal column (to flat sides)
set RcolDiag	[expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];	# Radius of octagonal column (to corners)
set Acol		[expr (2.0*$Dcol**2)/(1.0+sqrt(2.0))];			# Area of octagonal column section
set Jcol		[expr 1.2762*$RcolDiag**4];	    				# Polar moment of inertia for octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    				# Second moment of inertia, 1st transverse direction, for octagonal column section
set I3col		[expr 0.6381*$RcolDiag**4];    					# Second moment of inertia, 2nd transverse direction, for octagonal column section
set vecxzX 		[lindex $vecxzXList 12]; 						# X component of vecxz vector for local axis transformation
set vecxzY 		[lindex $vecxzYList 12]; 						# Y component of vecxz vector for local axis transformation
geomTransf $transformation 140 $vecxzX $vecxzY 0; 				# vecxz is perpendicular to the length of the cap beam.
if {$column_capbeam_joint == "rigidlink"} {
	set topN 1402
	set topS 1406
	set topM 1408
	rigidLink beam 1402 1403; 
	rigidLink beam 1406 1405; 
	rigidLink beam 1408 1404; 
} else {
	set topN 1403
	set topS 1405
	set topM 1404
}
if {$column_pins != 1 && $column_pins !=5} {
	make_column 14010 14010 $topN;
    make_column 14020 14070 $topS;
    make_column 14030 14090 $topM;
} elseif {$column_pins == 1} {
	make_column 14010 1401 $topN;
    make_column 14020 1407 $topS;
    make_column 14030 1409 $topM;
} elseif {$column_pins == 5} {
    set secTags 1390000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 14010;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	14010   1401   $topN    140  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
    set secTags 1390000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 14020;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	14020   1407   $topS    140  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
    set secTags 1390000
    for {set ip 1} {$ip < $np} {incr ip 1} {
        lappend secTags 14030;};
    set integration "FixedLocation $np $secTags $locations"
    element forceBeamColumn  	14030   1409   $topM    140  $integration; 	# Fixed Location integration points, point at BASE has PIN fiber section
}
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
set Acol		[expr (2*$Dcol**2)/(1+sqrt(2)) + ($Wcol-$Dcol)*$Dcol]; 	# Area of WIDE octagonal column section
set I2col		[expr 0.6381*$RcolDiag**4];	    						# Second moment of inertia about the 1st transverse (local y, horizontal) axis, for WIDE octagonal column section
set I3col		[expr $I2col+$Acol*$spO**2+($Wcol-$Dcol)*$Dcol**3/12]; 	# Second moment of inertia about the 2nd transverse (local z, vertical) axis, for WIDE octagonal column section
set Jcol		[expr $I2col+$I3col];		    						# Polar moment of inertia for WIDE octagonal column section
if {$column_capbeam_joint == "rigidlink"} {
	make_column 14040 1413 1414;
	rigidLink beam 1414 1415;
} else {
	make_column 14040 1413 1415;
}
# CAP BEAM ELEMENTS
# (Numbering scheme: 10000*Bent#+10,20,30,40)
# Transformation
# Deck and cap beam transformation, where vecxz is in the +Z direction
# (local z upwards on section, local y horiz left on section facing j (east))
geomTransf $transformation 301 0 0 1;
source ReadMPR.tcl; 										# Set up ReadMPR procedure for obtaining cap beam section properties
set CSDir "./Dimensions/CapCS/";  							# Directory containing cap beam cross section information
set CSType "Cap"; 											# Cross section type is cap beam
# Read in the concrete strength for each cap beam:
set fpfceCap [open "./Dimensions/fceCap.txt" r];
set fceCapList [read $fpfceCap];
close $fpfceCap
# Bents 2-10
for {set ib 2} {$ib <= 9} {incr ib 1} {
	lassign [ReadMPR $CSDir $CSType $ib {}] A Iy Iz J; # Cross Section properties
	set fceCap [lindex $fceCapList [expr $ib-2]]
	set EcCap [expr 57000.0*sqrt($fceCap*$ksi_psi)/$ksi_psi]
	set GcCap [expr $EcCap/(2.0*(1.0+$Uc))]
	element elasticBeamColumn [expr 10000*$ib+10] [expr 100*$ib+3] [expr 100*$ib+4] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	element elasticBeamColumn [expr 10000*$ib+20] [expr 100*$ib+4] [expr 100*$ib+5] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bents 10-11
for {set ib 10} {$ib <= 11} {incr ib 1} {
	lassign [ReadMPR $CSDir $CSType $ib {}] A Iy Iz J; # Cross Section properties
	set fceCap [lindex $fceCapList [expr $ib-2]]
	set EcCap [expr 57000.0*sqrt($fceCap*$ksi_psi)/$ksi_psi]
	set GcCap [expr $EcCap/(2.0*(1.0+$Uc))]
	element elasticBeamColumn [expr 10000*$ib+10] [expr 100*$ib+3] [expr 100*$ib+4] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	element elasticBeamColumn [expr 10000*$ib+20] [expr 100*$ib+4] [expr 100*$ib+5] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bent 12
lassign [ReadMPR $CSDir $CSType 12 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 120010 1203 1204 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120020 1204 1205 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120030 1205 1206 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120040 1206 1207 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 13 NE
lassign [ReadMPR $CSDir $CSType 13 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 130010 1303 1304 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 130020 1304 1305 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 14 NE
lassign [ReadMPR $CSDir $CSType 14 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 140010 1403 1404 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 140020 1404 1405 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# DECK ELEMENTS
# (Numbering scheme: 10*1stBent#+1stIntermediateNode#)
source ReadMPR.tcl; 										# Set up ReadMPR procedure for obtaining deck section properties
set CSDir "./Dimensions/DeckCS/";  							# Directory containing deck cross section information
set CSType "Deck"; 											# Cross section type is deck
# Abut 1 to Bent 2
lassign [ReadMPR $CSDir $CSType 1 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 10 1010 10001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 11 10001 10002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 12 10002 10003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 13 10003 10004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 14 10004 204 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 2 to Bent 5
for {set ib 2} {$ib <= 4} {incr ib 1} {
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 0}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib] [expr 100*$ib+4] [expr 10000*$ib+1] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 1}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+1] [expr 10000*$ib+1] [expr 10000*$ib+2] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 2}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+2] [expr 10000*$ib+2] [expr 10000*$ib+3] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 3}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+3] [expr 10000*$ib+3] [expr 10000*$ib+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 4}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+4] [expr 10000*$ib+4] [expr 100*($ib+1)+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bent 5 to Bent 6
lassign [ReadMPR $CSDir $CSType 5 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 50 504 50001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 5 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
if {$hinge_model != "none"} {
element elasticBeamColumn 51 500010 50002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
} else {
element elasticBeamColumn 51 50001 50002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
lassign [ReadMPR $CSDir $CSType 5 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 52 50002 50003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 5 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 53 50003 50004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 5 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 54 50004 604 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 6 to Bent 8
for {set ib 6} {$ib <= 7} {incr ib 1} {
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 0}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib] [expr 100*$ib+4] [expr 10000*$ib+1] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 1}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+1] [expr 10000*$ib+1] [expr 10000*$ib+2] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 2}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+2] [expr 10000*$ib+2] [expr 10000*$ib+3] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 3}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+3] [expr 10000*$ib+3] [expr 10000*$ib+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 4}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+4] [expr 10000*$ib+4] [expr 100*($ib+1)+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bent 8 to Bent 9
lassign [ReadMPR $CSDir $CSType 8 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 80 804 80001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 8 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 81 80001 80002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 8 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 82 80002 80003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 8 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 83 80003 80004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 8 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
if {$hinge_model != "none"} {
element elasticBeamColumn 84 800040 904 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
} else {
element elasticBeamColumn 84 80004 904 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;    
}
# Bent 9 to Bent 11
for {set ib 9} {$ib <= 10} {incr ib 1} {
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 0}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib] [expr 100*$ib+4] [expr 10000*$ib+1] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 1}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+1] [expr 10000*$ib+1] [expr 10000*$ib+2] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 2}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+2] [expr 10000*$ib+2] [expr 10000*$ib+3] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 3}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+3] [expr 10000*$ib+3] [expr 10000*$ib+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
	lassign [ReadMPR $CSDir $CSType $ib {-NodeNum 4}] A Iy Iz J; # Cross Section properties
	element elasticBeamColumn [expr 10*$ib+4] [expr 10000*$ib+4] [expr 100*($ib+1)+4] $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bent 11 to Bent 12
lassign [ReadMPR $CSDir $CSType 11 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 110 1104 110001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 11 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 111 110001 110002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 11 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 112 110002 110003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 11 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 113 110003 110004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 11 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 114 110004 1205 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 12 to Bent 13 NE
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 120 1204 120001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
if {$hinge_model != "none"} {
element elasticBeamColumn 121 1200010 120002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
} else {
element elasticBeamColumn 121 120001 120002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 122 120002 120003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 123 120003 120004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 124 120004 1304 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 12 to Bent 13 NR
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 5}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 125 1206 120005 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 6}] A Iy Iz J; # Cross Section properties
if {$hinge_model != "none"} {
element elasticBeamColumn 126 1200050 120006 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
} else {
element elasticBeamColumn 126 120005 120006 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 7}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 127 120006 120007 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 8}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 128 120007 120008 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 12 {-NodeNum 9}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 129 120008 1315 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 13 NE to Bent 14 NE
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 130 1304 130001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 131 130001 130002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 132 130002 130003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 133 130003 130004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 134 130004 1404 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 13 NR to Bent 14 NR
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 5}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 135 1315 130005 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 6}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 136 130005 130006 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 7}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 137 130006 130007 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 13 {-NodeNum 8}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 138 130007 130008 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 139 130008 1415 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 14 NE to Abut 15 NE
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 0}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 140 1404 140001 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 1}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 141 140001 140002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 2}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 142 140002 140003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 3}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 143 140003 140004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 4}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 144 140004 15010 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 14 NR to Abut 15 NR
lassign [ReadMPR $CSDir $CSType 14 {-NodeNum 5}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 145 1415 140005 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 146 140005 140006 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 147 140006 140007 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 148 140007 140008 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 149 140008 15040 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# ABUTMENT ELEMENTS
element elasticBeamColumn 1000 1020 1010 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-1
element elasticBeamColumn 1001 1010 1030 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-1
element elasticBeamColumn 15000 15020 15010 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-15NE
element elasticBeamColumn 15001 15010 15030 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-15NE
element elasticBeamColumn 15002 15050 15040 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-15NR
element elasticBeamColumn 15003 15040 15060 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;	# abut-15NR
#
# #+# COLUMN PINS ------------------------------------------------------------------
if {$column_capbeam_joint == "rigidlink"} {
	set topN 2
	set topS 6 
} else {
	set topN 3
	set topS 5
}
if {$column_pins == 2} {
	uniaxialMaterial Elastic 21		[expr 1.0e+15*$kips/$in]; 			# RIGID TRANSLATIONAL AND TORSIONAL STIFFNESS
	uniaxialMaterial Elastic 90 	[expr 0.0*$kips/$in]; 				# FREE X and/or Y ROTATIONAL STIFFNESS
	element zeroLength 	290000   201  2010  -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	290001   207  2070  -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	for {set ib 3} {$ib <= 11} {incr ib 1} {
		element zeroLength 	[expr $ib*100000+90000] [expr $ib*100+$topN] [expr $ib*1000+$topN*10] -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
		element zeroLength 	[expr $ib*100000+90001] [expr $ib*100+$topS] [expr $ib*1000+$topS*10] -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	}
	element zeroLength 	1290000  1201 12010 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1290001  1209 12090 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1290002  1211 12110 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1390000  1301 13010 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1390001  1307 13070 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1490000  1401 14010 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1490001  1407 14070 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
	element zeroLength 	1490002  1409 14090 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on;
}
if {$column_pins == 3} {
	uniaxialMaterial Elastic 21		[expr 1.0e+15*$kips/$in]; 			# RIGID TRANSLATIONAL AND TORSIONAL STIFFNESS
	if {$column_capbeam_joint == "none"} {
		uniaxialMaterial Elastic 90 	[expr 10.0*$kips/$in]; 			# FREE X and/or Y ROTATIONAL STIFFNESS
	} else {
		uniaxialMaterial Elastic 90 	[expr 0.0*$kips/$in]; 			# FREE X and/or Y ROTATIONAL STIFFNESS
	}
	element zeroLength 	290000   201  2010 									-mat 21 21 21 21 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 2 North Column, Base, DIR = {}
	element zeroLength 	290001   207  2070 									-mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 2 South Column, Base, DIR = {Y}
	element zeroLength 	390000   [expr 300+$topN]  	[expr 3000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 3 North Column, Top, DIR = {X}
	element zeroLength 	390001   [expr 300+$topS]  	[expr 3000+10*$topS] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 3 South Column, Top, DIR = {X,Y}
	element zeroLength 	490000   [expr 400+$topN]  	[expr 4000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 4 North Column, Top, DIR = {X}
	element zeroLength 	490001   [expr 400+$topS]  	[expr 4000+10*$topS] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 4 South Column, Top, DIR = {X,Y}
	element zeroLength 	590000   [expr 500+$topN]  	[expr 5000+10*$topN] 	-mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 5 North Column, Top, DIR = {Y}
	element zeroLength 	590001   [expr 500+$topS]  	[expr 5000+10*$topS] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 5 South Column, Top, DIR = {X}
	element zeroLength 	690000   [expr 600+$topN]  	[expr 6000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 6 North Column, Top, DIR = {X}
	element zeroLength 	690001   [expr 600+$topS]  	[expr 6000+10*$topS] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 6 South Column, Top, DIR = {X}
	element zeroLength 	790000   [expr 700+$topN]  	[expr 7000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 7 North Column, Top, DIR = {X}
	element zeroLength 	790001   [expr 700+$topS]  	[expr 7000+10*$topS] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 7 South Column, Top, DIR = {X,Y}
	element zeroLength 	890000   [expr 800+$topN]  	[expr 8000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 8 North Column, Top, DIR = {X}
	element zeroLength 	890001   [expr 800+$topS]  	[expr 8000+10*$topS] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 8 South Column, Top, DIR = {X,Y}
	element zeroLength 	990000   [expr 900+$topN]  	[expr 9000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 9 North Column, Top, DIR = {X}
	element zeroLength 	990001   [expr 900+$topS]  	[expr 9000+10*$topS] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 9 South Column, Top, DIR = {X}
	element zeroLength 	1090000  [expr 1000+$topN]  [expr 10000+10*$topN] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 10 North Column, Top, DIR = {X,Y}
	element zeroLength 	1090001  [expr 1000+$topS]  [expr 10000+10*$topS] 	-mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 10 South Column, Top, DIR = {X,Y}
	element zeroLength 	1190000  [expr 1100+$topN]  [expr 11000+10*$topN] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 11 North Column, Top, DIR = {X}
	element zeroLength 	1190001  [expr 1100+$topS]  [expr 11000+10*$topS] 	-mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 11 South Column, Top, DIR = {X}
	element zeroLength 	1290000  1201 12010 -mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 North Column, Base, DIR = {Y}
	element zeroLength 	1290001  1209 12090 -mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 South Column, Base, DIR = {Y}
	element zeroLength 	1290002  1211 12110 -mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 Center Column, Base, DIR = {Y}
	element zeroLength 	1390000  1301 13010 -mat 21 21 21 90 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 13 North Column, Base, DIR = {X}
	element zeroLength 	1390001  1307 13070 -mat 21 21 21 90 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 13 South Column, Base, DIR = {X,Y}
	element zeroLength 	1490000  1401 14010 -mat 21 21 21 21 21 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 North Column, Base, DIR = {}
	element zeroLength 	1490001  1407 14070 -mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 South Column, Base, DIR = {Y}
	element zeroLength 	1490002  1409 14090 -mat 21 21 21 21 90 21 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 Center Column, Base, DIR = {Y}
}
if {$column_pins == 4} {
	# Procedure for UNREINFORCED column fiber cross-section assignment (INELASTIC material properties):
	source ColSectionOctIEPin.tcl
	set Dcol [expr 84.0*$in];
	BuildOctColPINSection 290000  $Dcol; # Fiber cross-section for columns at Bents 2-11 (each have diameter 84.0 inches)
	set yp1 	[lindex $vecxzXList 0];
	set yp2 	[lindex $vecxzYList 0];
	element zeroLengthSection 290000 201 2010 290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 2 Left Column, Base, DIR = {}. Modeled as a zerolength octagonal section of unreinforced concrete.  Orientation vector sets the local x axis (perpendicular to the section) in the global Z (vertical) direction, and the local y axis parallel to the length of the cap beam.
	element zeroLengthSection 290001 207 2070 290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 2 Right Column, Base, DIR = {Y}
	for {set ib 3} {$ib <= 11} {incr ib 1} {
		set yp1 	[lindex $vecxzXList [expr $ib-2]];
		set yp2 	[lindex $vecxzYList [expr $ib-2]];
		element zeroLengthSection [expr $ib*100000+90000] [expr $ib*100+$topN] [expr $ib*1000+$topN*10] 290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent ib Left Column, Top
		element zeroLengthSection [expr $ib*100000+90001] [expr $ib*100+$topS] [expr $ib*1000+$topS*10] 290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent ib Right Column, Top
	}
	set Dcol [expr 66.0*$in];
	BuildOctColPINSection 1290000  $Dcol; # Fiber cross-section for columns at Bent 12 (each have diameter 66.0 inches) 
	set yp1 	[lindex $vecxzXList 10];
	set yp2 	[lindex $vecxzYList 10];
	element zeroLengthSection 1290000 1201 12010 1290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 Left Column, Base, DIR = {Y}
	element zeroLengthSection 1290001 1209 12090 1290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 Right Column, Base, DIR = {Y}
	element zeroLengthSection 1290002 1211 12110 1290000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 12 Center Column, Base, DIR = {Y}
	set Dcol [expr 48.0*$in];
	BuildOctColPINSection 1390000  $Dcol; # Fiber cross-section for columns at Bents 13-14 (each have diameter 48.0 inches)
	set yp1 	[lindex $vecxzXList 11];
	set yp2 	[lindex $vecxzYList 11];
	element zeroLengthSection 1390000 1301 13010 1390000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 13 Left Column, Base, DIR = {X}
	element zeroLengthSection 1390001 1307 13070 1390000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 13 Right Column, Base, DIR = {X,Y}
	set yp1 	[lindex $vecxzXList 12];
	set yp2 	[lindex $vecxzYList 12];
	element zeroLengthSection 1490000 1401 14010 1390000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 Left Column, Base, DIR = {}
	element zeroLengthSection 1490001 1407 14070 1390000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 Right Column, Base, DIR = {Y}
	element zeroLengthSection 1490002 1409 14090 1390000 -orient 0 0 1 $yp1 $yp2 0 -doRayleigh $rayleigh_zerolength_on; # PIN Bent 14 Center Column, Base, DIR = {Y}
}
#
# #+# IN-SPAN HINGE SPRINGS --------------------------------------------------------
if {$hinge_model != "none"} {
set khlc 						[expr 25000.0*$kips/$in];  			# Pounding stiffness is equal for all hinges' longitudinal compressive stiffness
set phl 						[expr 1.0e+9*$kips];     			# Pounding strength and tensile strength are both a high number for all hinges' longitudinal strength
set Eb 							[expr 0.50*$ksi]; 					# Elastic modulus for elastomeric bearings
set Gb 							[expr $CGh*0.1*$ksi]; 				# Shear modulus for elastomeric bearings
set gapV	    				[expr 0.6*$in];  					# The flexible portion of the elastomeric bearing pad
# In span hinge 5 Properties
set lb5 						[expr 4.5*$in]; 					# Hinge 5 elastomeric bearing thickness
set ab5 						[expr 16.0*22.0*$in**2];			# Hinge 5 elastomeric bearing area
set nb5 						[expr 6.0]; 						# Hinge 5 number of elastomeric bearings
set khlb5 						[expr $nb5*$Gb*$ab5/$lb5]; 			# Hinge 5 longitudinal initial stiffness from shear resistance of bearing pads (in kip/in) 
set gaphl5 						[expr 3.50*$in]; 					# Longitudinal gap width Hinge 5  
set alr5 						[expr 24.54*$in**2]; 				# Longitudinal restrainer area Hinge 5
set llr5 						[expr 506.0*$in]; 					# Longitudinal restrainer length Hinge 5
set khlt5 						[expr $Es*$alr5/$llr5]; 			# Longitudinal tensile stiffness Hinge 5
set kht5 						[expr 41422.0*$kips/$in]; 			# Transverse stiffness Hinge 5
set pht5 						[expr 1036.0*$kips]; 				# Transverse strength Hinge 5
set khv5 						[expr $nb5*$Eb*$ab5/$lb5]; 			# Hinge 5 vertical initial stiffness from compressive resistance of bearing pads (in kip/in)
# In span hinge 8 Properties
set lb8 						[expr 3.5*$in]; 					# Hinge 8 elastomeric bearing thickness
set ab8 						[expr 12.0*22.0*$in**2];			# Hinge 8 elastomeric bearing area
set nb8 						[expr 6.0]; 						# Hinge 8 number of elastomeric bearings
set khlb8 						[expr $nb8*$Gb*$ab8/$lb8]; 			# Hinge 8 longitudinal initial stiffness from shear resistance of bearing pads (in kip/in) 
set gaphl8 						[expr 2.75*$in]; 					# Longitudinal gap width Hinge 8
set alr8 						[expr 24.54*$in**2]; 				# Longitudinal restrainer area Hinge 8
set llr8 						[expr 458.0*$in]; 					# Longitudinal restrainer length Hinge 8
set khlt8 						[expr $Es*$alr8/$llr8]; 			# Longitudinal tensile stiffness Hinge 8
set kht8 						[expr 41422.0*$kips/$in]; 			# Transverse stiffness Hinge 8
set pht8 						[expr 1036.0*$kips]; 				# Transverse strength Hinge 8
set khv8 						[expr $nb8*$Eb*$ab8/$lb8]; 			# Hinge 8 vertical initial stiffness from compressive resistance of bearing pads (in kip/in)
# In span hinge 12 NE Properties
set lb12NE 						[expr 3.0*$in]; 					# Hinge 12 NE elastomeric bearing thickness
set ab12NE 						[expr 14.0*22.0*$in**2];			# Hinge 12 NE elastomeric bearing area
set nb12NE 						[expr 6.0]; 						# Hinge 12 NE number of elastomeric bearings
set khlb12NE 					[expr $nb12NE*$Gb*$ab12NE/$lb12NE]; # Hinge 12 NE longitudinal initial stiffness from shear resistance of bearing pads (in kip/in)
set gaphl12NE 					[expr 2.25*$in]; 					# Longitudinal gap width Hinge 12 NE
set alr12NE 					[expr 24.54*$in**2]; 				# Longitudinal restrainer area Hinge 12 NE
set llr12NE 					[expr 326.0*$in]; 					# Longitudinal restrainer length Hinge 12 NE
set khlt12NE 					[expr $Es*$alr12NE/$llr12NE]; 		# Longitudinal tensile stiffness Hinge 12 NE
set kht12NE 					[expr 41422.0*$kips/$in]; 			# Transverse stiffness Hinge 12 NE
set pht12NE 					[expr 1036.0*$kips]; 				# Transverse strength Hinge 12 NE
set khv12NE 					[expr $nb12NE*$Eb*$ab12NE/$lb12NE]; # Hinge 12 NE vertical initial stiffness from compressive resistance of bearing pads (in kip/in)
# In span hinge 12 NR Properties
set lb12NR 						[expr 2.5*$in]; 					# Hinge 12 NR elastomeric bearing thickness
set ab12NR 						[expr 12.0*18.0*$in**2];			# Hinge 12 NR elastomeric bearing area
set nb12NR 						[expr 3.0]; 						# Hinge 12 NR number of elastomeric bearings
set khlb12NR 					[expr $nb12NR*$Gb*$ab12NR/$lb12NR]; # Hinge 12 NR longitudinal initial stiffness from shear resistance of bearing pads (in kip/in)
set gaphl12NR 					[expr 1.75*$in]; 					# Longitudinal gap width Hinge 12 NR
set alr12NR 					[expr 9.82*$in**2]; 				# Longitudinal restrainer area Hinge 12 NR
set llr12NR 					[expr 330.0*$in]; 					# Longitudinal restrainer length Hinge 12 NR
set khlt12NR 					[expr $Es*$alr12NR/$llr12NR]; 		# Longitudinal tensile stiffness Hinge 12 NR
set kht12NR 					[expr 34518.0*$kips/$in]; 			# Transverse stiffness Hinge 12 NR
set pht12NR 					[expr 863.0*$kips]; 				# Transverse strength Hinge 12 NR
set lb12NR 						[expr 2.5*$in]; 					# Hinge 12 NR elastomeric bearing thickness
set ab12NR 						[expr 12.0*18.0*$in**2];			# Hinge 12 NR elastomeric bearing area
set nb12NR 						[expr 3.0]; 						# Hinge 12 NR number of elastomeric bearings
set khv12NR 					[expr $nb12NR*$Eb*$ab12NR/$lb12NR]; # Hinge 12 NR vertical initial stiffness from compressive resistance of bearing pads (in kip/in)
# puts "\tkhlc = $khlc; phl = $phl"
# puts "\tkhlt5 = $khlt5; khlt8 = $khlt8; khlt12NE = $khlt12NE; khlt12NR = $khlt12NR"
# puts "\tkhlb5 = $khlb5; khlb8 = $khlb8; khlb12NE = $khlb12NE; khlb12NR = $khlb12NR"
# puts "\tkht5 = $kht5; kht8 = $kht8; kht12NE = $kht12NE, kht12NR = $kht12NR"
# puts "\tpht5 = $pht5; pht8 = $pht8; pht12NE = $pht12NE, pht12NR = $pht12NR"
# puts "\tkhv5 = $khv5; khv8 = $khv8; khv12NE = $khv12NE, khv12NR = $khv12NR"
# General release (zero stiffness) spring for rotational dofs
uniaxialMaterial Elastic 		504 	[expr 0.0*$kips/$in]; 																# ROTATIONAL stiffness for all hinges set to zero
if {$hinge_model == "linear"} {
# Hinge 5 Spring Materials
uniaxialMaterial Elastic 		501 	$khlb5; 																			# LONGITUDINAL Hinge 5 (BEARING SHEAR)
uniaxialMaterial Elastic 		502 	$kht5; 																				# TRANSVERSE Hinge 5 (SHEAR KEY)
uniaxialMaterial Elastic 		503 	$khv5; 																				# VERTICAL Hinge 5 (BEARING COMPRESSION)
# Hinge 8 Spring Materials
uniaxialMaterial Elastic 		801 	$khlb8; 																			# LONGITUDINAL Hinge 8 (BEARING SHEAR)
uniaxialMaterial Elastic 		802 	$kht8; 																				# TRANSVERSE Hinge 8 (SHEAR KEY)
uniaxialMaterial Elastic 		803 	$khv8; 																				# VERTICAL Hinge 8 (BEARING COMPRESSION)
# Hinge 12 NE Spring Materials
uniaxialMaterial Elastic 		1201 	$khlb12NE; 																			# LONGITUDINAL Hinge 12 NE (BEARING SHEAR)
uniaxialMaterial Elastic 		1202 	$kht12NE; 																			# TRANSVERSE Hinge 12 NE (SHEAR KEY)
uniaxialMaterial Elastic 		1203 	$khv12NE; 																			# VERTICAL Hinge 12 NE (BEARING COMPRESSION)
# Hinge 12 NR Spring Materials
uniaxialMaterial Elastic 		1204 	$khlb12NR; 																			# LONGITUDINAL Hinge 12 NR (BEARING SHEAR)
uniaxialMaterial Elastic 		1205 	$kht12NR; 																			# TRANSVERSE Hinge 12 NR (SHEAR KEY)
uniaxialMaterial Elastic 		1206 	$khv12NR; 																			# VERTICAL Hinge 12 NR (BEARING COMPRESSION)
}
if {$hinge_model == "simplified"} {
# Hinge 5 Spring Materials
uniaxialMaterial ElasticPPGap 	5011 	[expr $khlc] [expr -$phl] [expr -$gaphl5] 1.0e-3 damage; 							# Hinge 5 longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	5012 	[expr $khlt5] [expr $phl] [expr $gaphl5] 1.0e-3 damage; 							# Hinge 5 longitudinal TENSION
uniaxialMaterial Elastic 		5013 	$khlb5; 																			# Hinge 5 longitudinal BEARING SHEAR
uniaxialMaterial Parallel 		501 	5011 5012 5013; 																	# LONGITUDINAL Hinge 5
uniaxialMaterial Elastic 		502 	$kht5; 																				# TRANSVERSE Hinge 5 (SHEAR KEY)
uniaxialMaterial ENT 			503 	$khv5; 																				# VERTICAL Hinge 5 (BEARING COMPRESSION)
# Hinge 8 Spring Materials
uniaxialMaterial ElasticPPGap 	8011 	[expr $khlc] [expr -$phl] [expr -$gaphl8] 1.0e-3 damage; 							# Hinge 8 longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	8012 	[expr $khlt8] [expr $phl] [expr $gaphl8] 1.0e-3 damage; 							# Hinge 8 longitudinal TENSION
uniaxialMaterial Elastic 		8013 	$khlb8; 																			# Hinge 8 longitudinal BEARING SHEAR
uniaxialMaterial Parallel 		801 	8011 8012 8013; 																	# LONGITUDINAL Hinge 8
uniaxialMaterial Elastic 		802 	$kht8; 																				# TRANSVERSE Hinge 8 (SHEAR KEY)
uniaxialMaterial ENT 			803 	$khv8; 																				# VERTICAL Hinge 8 (BEARING COMPRESSION)
# Hinge 12 NE Spring Materials
uniaxialMaterial ElasticPPGap 	12011 	[expr $khlc] [expr -$phl] [expr -$gaphl12NE] 1.0e-3 damage; 						# Hinge 12 NE longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	12012 	[expr $khlt12NE] [expr $phl] [expr $gaphl12NE] 1.0e-3 damage; 						# Hinge 12 NE longitudinal TENSION
uniaxialMaterial Elastic 		12013 	$khlb12NE; 																			# Hinge 12 NE longitudinal BEARING SHEAR
uniaxialMaterial Parallel 		1201 	12011 12012 12013; 																	# LONGITUDINAL Hinge 12 NE
uniaxialMaterial Elastic 		1202 	$kht12NE; 																			# TRANSVERSE Hinge 12 NE (SHEAR KEY)
uniaxialMaterial ENT 			1203 	$khv12NE; 																			# VERTICAL Hinge 12 NE (BEARING COMPRESSION)
# Hinge 12 NR Spring Materials
uniaxialMaterial ElasticPPGap 	12041 	[expr $khlc] [expr -$phl] [expr -$gaphl12NR] 1.0e-3 damage; 						# Hinge 12 NR longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	12042 	[expr $khlt12NR] [expr $phl] [expr $gaphl12NR] 1.0e-3 damage; 						# Hinge 12 NR longitudinal TENSION
uniaxialMaterial Elastic 		12043 	$khlb12NR; 																			# Hinge 12 NR longitudinal BEARING SHEAR
uniaxialMaterial Parallel 		1204 	12041 12042 12043; 																	# LONGITUDINAL Hinge 12 NR
uniaxialMaterial Elastic 		1205 	$kht12NR; 																			# TRANSVERSE Hinge 12 NR (SHEAR KEY)
uniaxialMaterial ENT 			1206 	$khv12NR; 																			# VERTICAL Hinge 12 NR (BEARING COMPRESSION)
}
if {$hinge_model == "complex"} {
# Hinge 5 Spring Materials
uniaxialMaterial ElasticPPGap 	5011 	[expr $khlc] [expr -$phl] [expr -$gaphl5] 1.0e-3 damage; 							# Hinge 5 longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	5012 	[expr $khlt5] [expr $phl] [expr $gaphl5] 1.0e-3 damage; 							# Hinge 5 longitudinal TENSION
uniaxialMaterial Elastic 		5013 	$khlb5; 																			# Hinge 5 longitudinal BEARING SHEAR
uniaxialMaterial Parallel 		501 	5011 5012 5013; 																	# LONGITUDINAL Hinge 5
uniaxialMaterial ElasticPP 		502 	$kht5 [expr $pht5/$kht5]; 															# TRANSVERSE Hinge 5 (SHEAR KEY)
uniaxialMaterial ENT 			5031 	$khv5; 																				# Hinge 5 vertical bearing compression
uniaxialMaterial ElasticPPGap 	5032 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage; 								# Hinge 5 vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		503 	5031 5032; 																			# VERTICAL Hinge 5
# Hinge 8 Spring Materials
uniaxialMaterial ElasticPPGap 	8011 	[expr $khlc] [expr -$phl] [expr -$gaphl8] 1.0e-3 damage; 							# Hinge 8 longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	8012 	[expr $khlt8] [expr $phl] [expr $gaphl8] 1.0e-3 damage; 							# Hinge 8 longitudinal TENSION
uniaxialMaterial Elastic 		8013 	$khlb8; 																			# Hinge 8 BEARING SHEAR
uniaxialMaterial Parallel 		801 	8011 8012 8013; 																	# LONGITUDINAL Hinge 8
uniaxialMaterial ElasticPP 		802 	$kht8 [expr $pht8/$kht8]; 															# TRANSVERSE Hinge 8 (SHEAR KEY)
uniaxialMaterial ENT 			8031 	$khv8; 																				# Hinge 8 vertical bearing compression
uniaxialMaterial ElasticPPGap 	8032 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage; 								# Hinge 8 vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		803 	8031 8032; 																			# VERTICAL Hinge 8
# Hinge 12 NE Spring Materials
uniaxialMaterial ElasticPPGap 	12011 	[expr $khlc] [expr -$phl] [expr -$gaphl12NE] 1.0e-3 damage; 						# Hinge 12 NE longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	12012 	[expr $khlt12NE] [expr $phl] [expr $gaphl12NE] 1.0e-3 damage; 						# Hinge 12 NE longitudinal TENSION
uniaxialMaterial Elastic 		12013 	$khlb12NE; 																			# Hinge 12 NE BEARING SHEAR
uniaxialMaterial Parallel 		1201 	12011 12012 12013; 																	# LONGITUDINAL Hinge 12 NE
uniaxialMaterial ElasticPP 		1202 	$kht12NE [expr $pht12NE/$kht12NE]; 													# TRANSVERSE Hinge 12 NE (SHEAR KEY)
uniaxialMaterial ENT 			12031 	$khv12NE; 																			# Hinge 12 NE vertical bearing compression
uniaxialMaterial ElasticPPGap 	12032 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage;									# Hinge 12 NE vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		1203 	12031 12032; 																		# VERTICAL Hinge 12 NE
# Hinge 12 NR Spring Materials
uniaxialMaterial ElasticPPGap 	12041 	[expr $khlc] [expr -$phl] [expr -$gaphl12NR] 1.0e-3 damage; 						# Hinge 12 NR longitudinal COMPRESSION
uniaxialMaterial ElasticPPGap 	12042 	[expr $khlt12NR] [expr $phl] [expr $gaphl12NR] 1.0e-3 damage; 						# Hinge 12 NR longitudinal TENSION
uniaxialMaterial Elastic 		12043 	$khlb12NR; 																			# Hinge 12 NR BEARING SHEAR
uniaxialMaterial Parallel 		1204 	12041 12042 12043; 																	# LONGITUDINAL Hinge 12 NR
uniaxialMaterial ElasticPP 		1205 	$kht12NR [expr $pht12NR/$kht12NR]; 													# TRANSVERSE Hinge 12 NR (SHEAR KEY)
uniaxialMaterial ENT 			12061 	$khv12NR; 																			# Hinge 12 NR vertical bearing compression
uniaxialMaterial ElasticPPGap 	12062 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage;									# Hinge 12 NR vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		1206 	12061 12062; 																		# VERTICAL Hinge 12 NR
}
# ZeroLength elements defined for each hinge
element zeroLength 				500000  50001 500010 -mat 501 502 503 504 504 504 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 456.63 -162.15534 0 0 1 0; # IN-SPAN HINGE BENT 5
element zeroLength 				800000  80004 800040 -mat 801 802 803 504 504 504 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 702.52 -66.80943 0 0 -1 0; # IN-SPAN HINGE BENT 8 (Needs -Y as the yp vector because the second node is the support)
element zeroLength 				1200000 120001 1200010 -mat 1201 1202 1203 504 504 504 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 353.11 41.662 0 0 1 0; # IN-SPAN HINGE BENT 12 NE
element zeroLength 				1200001 120005 1200050 -mat 1204 1205 1206 504 504 504 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 195.72 -6.42 0 0 1 0; # IN-SPAN HINGE BENT 12 NR
}
#
# #+# ABUTMENT SPRINGS -------------------------------------------------------------
if {$abutment_model != "none"} {
set CW 			[expr 4.0/3.0]; 		# Wall participation coefficient (for wingwall contribution to transverse stiffness)
set CL 			[expr 2.0/3.0]; 		# Wall effectiveness (for wingwall contribution to transverse stiffness)
set CWT 		[expr 1.0/3.0];  		# Multiplier for assumed effective wingwall width = 1/3 of backwall width
set Eb 			[expr 0.50*$ksi]; 		# Elastic modulus for elastomeric bearings
set Gb 			[expr $CGa*0.1*$ksi]; 	# Shear modulus for elastomeric bearings
set gapV	    [expr 0.6*$in];  		# The plastic portion of the elastomeric bearing pad which is flexible
# Abutment 1 Properties
set wabut1 		[expr 54.0]; 											# Abutment 1 backwall width (in feet)
set habut1 		[expr 7.0]; 											# Abutment 1 backwall height (in feet)
set theta1 		[expr 0.25]; 											# Abutment 1 skew angle (in degrees)
set Rsk1 		[expr exp(-$theta1/45)]; 								# Abutment 1 skew reduction factor
set kabutl1 	[expr $wabut1*(5.5*$habut1+20.0)*$Rsk1*$kips/$in]; 		# Abutment 1 longitudinal stiffness (in kip/in)
set kabutt1 	[expr $CW*$CL*$CWT*$kabutl1*$kips/$in]; 				# Abutment 1 transverse stiffness (in kip/in)
set pabutl1 	[expr $wabut1*(5.5*$habut1**2.5)*$Rsk1/(1+2.37*$habut1)];# Abutment 1 longitudinal resistance force (in kips)
set pabutt1 	[expr $CW*$CL*$CWT*$pabutl1*$kips]; 					# Abutment 1 transverse resistance force (in kips)
set gap1 		[expr 1.75*$in];  										# Abutment 1 longitudinal gap (in inches)
set lb1 		[expr 2.0*$in]; 										# Abutment 1 elastomeric bearing thickness
set ab1 		[expr 12.0*18.0*$in**2];								# Abutment 1 elastomeric bearing area
set nb1 		[expr 6.0]; 											# Abutment 1 number of elastomeric bearings
set kabutlb1 	[expr $nb1*$Gb*$ab1/$lb1]; 			                    # Abutment 1 longitudinal initial stiffness from shear resistance of bearing pads (in kip/in) 
set kabutv1 	[expr $nb1*$Eb*$ab1/$lb1]; 								# Abutment 1 vertical initial stiffness (in kip/in)
# Abutment 15 NE Properties
set wabut15NE 	[expr 75.5]; 											# Abutment 15 NE width (in feet)
set habut15NE 	[expr 5.5]; 											# Abutment 15 NE backwall height (in feet)
set theta15NE 	[expr 44.34]; 											# Abutment 15 NE skew angle (in degrees)
set Rsk15NE 	[expr exp(-$theta15NE/45)]; 							# Abutment 15 NE skew reduction factor
set kabutl15NE 	[expr $wabut15NE*(5.5*$habut15NE+20.0)*$Rsk15NE*$kips/$in];# Abutment 15 NE longitudinal stiffness (in kip/in)
set kabutt15NE 	[expr $CW*$CL*$CWT*$kabutl15NE*$kips/$in]; 				# Abutment 15 NE transverse stiffness (in kip/in)
set pabutl15NE 	[expr $wabut15NE*(5.5*$habut15NE**2.5)*$Rsk15NE/(1+2.37*$habut15NE)];# Abutment 15 NE longitudinal resistance force (in kips)
set pabutt15NE 	[expr $CW*$CL*$CWT*$pabutl15NE*$kips]; 					# Abutment 15 NE transverse resistance force (in kips)
set gap15NE 	[expr 1.0*$in];  										# Abutment 15 NE longitudinal gap (in inches)
set lb15NE 		[expr 1.50*$in]; 										# Abutment 15 NE elastomeric bearing thickness
set ab15NE 		[expr 16.0*20.0*$in**2];								# Abutment 15 NE elastomeric bearing area
set nb15NE 		[expr 6.0]; 											# Abutment 15 NE number of elastomeric bearings
set kabutlb15NE [expr $nb15NE*$Gb*$ab15NE/$lb15NE]; 			        # Abutment 15 NE longitudinal initial stiffness from shear resistance of bearing pads (in kip/in) 
set kabutv15NE 	[expr $nb15NE*$Eb*$ab15NE/$lb15NE]; 					# Abutment 15 NE vertical initial stiffness (in kip/in)
# Abutment 15 NR Properties
set wabut15NR 	[expr 32.1]; 											# Abutment 15 NR width (in feet)
set habut15NR 	[expr 5.5]; 											# Abutment 15 NR backwall height (in feet)
set theta15NR 	[expr 35.896]; 											# Abutment 15 NR skew angle (in degrees)
set Rsk15NR 	[expr exp(-$theta15NR/45)]; 							# Abutment 15 NR skew reduction factor
set kabutl15NR 	[expr $wabut15NR*(5.5*$habut15NR+20.0)*$Rsk15NR]; 		# Abutment 15 NR longitudinal stiffness (in kip/in)
set kabutt15NR 	[expr $CW*$CL*$CWT*$kabutl15NR*$kips/$in]; 				# Abutment 15 NR transverse stiffness (in kip/in)
set pabutl15NR 	[expr $wabut15NR*(5.5*$habut15NR**2.5)*$Rsk15NR/(1+2.37*$habut15NR)];# Abutment 15 NR longitudinal resistance force (in kips)
set pabutt15NR 	[expr $CW*$CL*$CWT*$pabutl15NR*$kips]; 					# Abutment 15 NR transverse resistance force (in kips)
set gap15NR 	[expr 1.0*$in];  										# Abutment 15 NR longitudinal gap (in inches)
set lb15NR 		[expr 1.50*$in]; 										# Abutment 15 NR elastomeric bearing thickness
set ab15NR 		[expr 16.0*20.0*$in**2];								# Abutment 15 NR elastomeric bearing area
set nb15NR 		[expr 3.0]; 											# Abutment 15 NR number of elastomeric bearings
set kabutlb15NR [expr $nb15NR*$Gb*$ab15NR/$lb15NR]; 			        # Abutment 15 NR longitudinal initial stiffness from shear resistance of bearing pads (in kip/in) 
set kabutv15NR 	[expr $nb15NR*$Eb*$ab15NR/$lb15NR]; 					# Abutment 15 NR vertical initial stiffness (in kip/in)
# puts "\tkabutlb1 = $kabutlb1; kabutlb15NE = $kabutlb15NE; kabutlb15NR = $kabutlb15NR"
# puts "\tkabutl1 = $kabutl1; kabutl15NE = $kabutl15NE; kabutl15NR = $kabutl15NR"
# puts "\tkabutt1 = $kabutt1; kabutt15NE = $kabutt15NE; kabutt15NR = $kabutt15NR"
# puts "\tkabutv1 = $kabutv1; kabutv15NE = $kabutv15NE; kabutv15NR = $kabutv15NR"
# puts "\tpabutl1 = $pabutl1; pabutl15NE = $pabutl15NE, pabutl15NR = $pabutl15NR"
# puts "\tpabutt1 = $pabutt1, pabutt15NE = $pabutt15NE, pabutt15NR = $pabutt15NR"
set cfactor	 	  0.5;
# Spring Material for all abutments
uniaxialMaterial Elastic 		2000 	[expr 0.0*$kips/$in]; 																		# ROTATIONAL stiffness for all abutments set to zero
if {$abutment_model == "linear"} {		
# Abutment 1 Spring Materials		
uniaxialMaterial Elastic 		201	    $kabutlb1;  																				# LONGITUDINAL Abutment 1 (BEARING SHEAR)
uniaxialMaterial Elastic 		202 	[expr $kabutt1*$cfactor]; 								                                    # TRANSVERSE Abutment 1 (SHEAR KEY)
uniaxialMaterial Elastic 		203 	$kabutv1;																					# VERTICAL Abutment 1 (BEARING COMPRESSION)
# Abutment 15 NE Spring Materials	
uniaxialMaterial Elastic 		215	    $kabutlb15NE;  																				# LONGITUDINAL Abutment 15 NE (BEARING SHEAR)
uniaxialMaterial Elastic 		216 	[expr $kabutt15NE*$cfactor]; 								                                # TRANSVERSE Abutment 15 NE (SHEAR KEY)
uniaxialMaterial Elastic 		217 	$kabutv15NE;																				# VERTICAL Abutment 15 NE (BEARING COMPRESSION)
# Abutment 15 NR Spring Materials		
uniaxialMaterial Elastic 		218	    $kabutlb15NR;																				# LONGITUDINAL Abutment 15 NR (BEARING SHEAR)
uniaxialMaterial Elastic 		219 	[expr $kabutt15NR*$cfactor]; 								                                # TRANSVERSE Abutment 15 NR (SHEAR KEY)
uniaxialMaterial Elastic 		220 	$kabutv15NR;																				# VERTICAL Abutment 15 NR (BEARING COMPRESSION)
}
if {$abutment_model == "simplified"} {
# Abutment 1 Spring Materials
uniaxialMaterial ElasticPPGap 	2011 	[expr $kabutl1*$cfactor] [expr -$pabutl1*$cfactor] [expr -$gap1] 1.0e-3 damage; 			# Abutment 1 wingwall resistance
uniaxialMaterial Elastic 		2012 	$kabutlb1; 																					# Abutment 1 bearing shear
uniaxialMaterial Parallel       201	    2011 2012;  																				# LONGITUDINAL Abutment 1
uniaxialMaterial Elastic 		202 	[expr $kabutt1*$cfactor]; 								                                    # TRANSVERSE Abutment 1 (SHEAR KEY)
uniaxialMaterial ENT 			203 	$kabutv1;																					# VERTICAL Abutment 1 (BEARING COMPRESSION)
# Abutment 15 NE Spring Materials
uniaxialMaterial ElasticPPGap 	2151 	[expr $kabutl15NE*$cfactor] [expr -$pabutl15NE*$cfactor] [expr -$gap15NE] 1.0e-3 damage; 	# Abutment 15 NE wingwall resistance
uniaxialMaterial Elastic 		2152 	$kabutlb15NE; 																				# Abutment 15 NE bearing shear
uniaxialMaterial Parallel       215	    2151 2152;  																				# LONGITUDINAL Abutment 15 NE
uniaxialMaterial Elastic 		216 	[expr $kabutt15NE*$cfactor]; 								                                # TRANSVERSE Abutment 15 NE (SHEAR KEY)
uniaxialMaterial ENT 			217 	$kabutv15NE;																				# VERTICAL Abutment 15 NE (BEARING COMPRESSION)
# Abutment 15 NR Spring Materials
uniaxialMaterial ElasticPPGap 	2181 	[expr $kabutl15NR*$cfactor] [expr -$pabutl15NR*$cfactor] [expr -$gap15NR] 1.0e-3 damage; 	# Abutment 15 NR wingwall resistance
uniaxialMaterial Elastic 		2182 	$kabutlb15NR; 																				# Abutment 15 NR bearing shear
uniaxialMaterial Parallel       218	    2181 2182; 																					# LONGITUDINAL Abutment 15 NR
uniaxialMaterial Elastic 		219 	[expr $kabutt15NR*$cfactor]; 								                                # TRANSVERSE Abutment 15 NR (SHEAR KEY)
uniaxialMaterial ENT 			220 	$kabutv15NR;																				# VERTICAL Abutment 15 NR (BEARING COMPRESSION)
}
if {$abutment_model == "complex"} {
# Abutment 1 Spring Materials
uniaxialMaterial ElasticPPGap 	2011 	[expr $kabutl1*$cfactor] [expr -$pabutl1*$cfactor] [expr -$gap1] 1.0e-3 damage; 			# Abutment 1 wingwall resistance
uniaxialMaterial Elastic 		2012 	$kabutlb1; 																					# Abutment 1 bearing shear
uniaxialMaterial Parallel       201	    2011 2012;  																				# LONGITUDINAL Abutment 1
uniaxialMaterial ElasticPP 		202 	[expr $kabutt1*$cfactor] [expr $pabutt1/$kabutt1]; 											# TRANSVERSE Abutment 1
uniaxialMaterial ENT 			2031 	$kabutv1; 																				# Abutment 1 vertical bearing compression
uniaxialMaterial ElasticPPGap 	2032 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage; 										# Abutment 1 vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		203 	2031 2032; 																					# VERTICAL Abutment 1
# Abutment 15 NE Spring Materials
uniaxialMaterial ElasticPPGap 	2151 	[expr $kabutl15NE*$cfactor] [expr -$pabutl15NE*$cfactor] [expr -$gap15NE] 1.0e-3 damage; 	# Abutment 15 NE wingwall resistance
uniaxialMaterial Elastic 		2152 	$kabutlb15NE; 																				# Abutment 15 NE bearing shear
uniaxialMaterial Parallel       215	    2151 2152;  																				# LONGITUDINAL Abutment 15 NE
uniaxialMaterial ElasticPP 		216 	[expr $kabutt15NE*$cfactor] [expr $pabutt15NE/$kabutt15NE]; 								# TRANSVERSE Abutment 15 NE
uniaxialMaterial ENT 			2171 	$kabutv15NE; 																				# Abutment 15 NE vertical bearing compression
uniaxialMaterial ElasticPPGap 	2172 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage; 										# Abutment 15 NE vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		217 	2171 2172; 																					# VERTICAL Abutment 15 NE
# Abutment 15 NR Spring Materials
uniaxialMaterial ElasticPPGap 	2181 	[expr $kabutl15NR*$cfactor] [expr -$pabutl15NR*$cfactor] [expr -$gap15NR] 1.0e-3 damage; 	# Abutment 15 NR wingwall resistance
uniaxialMaterial Elastic 		2182 	$kabutlb15NR; 																				# Abutment 15 NR bearing shear
uniaxialMaterial Parallel       218	    2181 2182; 																					# LONGITUDINAL Abutment 15 NR
uniaxialMaterial ElasticPP 		219 	[expr $kabutt15NR*$cfactor] [expr $pabutt15NR/$kabutt15NR]; 								# TRANSVERSE Abutment 15 NR
uniaxialMaterial ENT 			2201 	$kabutv15NR; 																				# Abutment 15 NR vertical bearing compression
uniaxialMaterial ElasticPPGap 	2202 	[expr 1.0e+9] [expr -1.0e+10] -$gapV 1.0e-3 damage; 										# Abutment 15 NR vertical rigid zone after bearing is compressed
uniaxialMaterial Parallel 		220 	2201 2202; 																					# VERTICAL Abutment 15 NR
}
# zerolength element defined from fixed end to the free end of the spring and the z axis defined by cross product of "x" and "yp"
# to have the rigid motion of the abutment, the rotational degree of freedom of the zerolength element is released
# Note that the local x axis is parallel to the backwall (abutment element), which corresponds to the transverse direction. The local y axis is perpendicular to the backwall and corresponds to the longitudinal direction.
element zeroLength 120000  1021 1020 -mat 202 201 203 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient -156.24 -202.12 0 176.95 -135.55 0; 
element zeroLength 130000  1031 1030 -mat 202 201 203 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 156.24 202.12 0 -176.95 135.55 0; 
element zeroLength 1510000  15021 15020 -mat 216 215 217 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient -223.99 -290.74 0 329.51 38.87 0; 
element zeroLength 1520000  15031 15030 -mat 216 215 217 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 223.99 290.74 0 -329.51 -38.87 0; 
element zeroLength 1550000  15051 15050 -mat 219 218 220 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient -77.97 -101.20 0 348.96 -10.50 0; 
element zeroLength 1560000  15061 15060 -mat 219 218 220 2000 2000 2000 -dir 1 2 3 4 5 6 -doRayleigh $rayleigh_zerolength_on -orient 77.97 101.20 0 -348.96 10.50 0;
}
#
# #+# TIMER 1 ----------------------------------------------------------------------
# Record time post model definition
set tModel [expr {([clock clicks -millisec]-$t0)/1000.}];
puts stderr "\tThe model definition time was $tModel seconds"; # Print out the model definition time
puts $timers "{\"tModel\": $tModel,";
set t1 [clock clicks -millisec];  
#
# #+# GRAVITY RECORDERS ------------------------------------------------------------
set gravRec ""
# DISPLACEMENTS AT ALL NODES
lappend gravRec [recorder Node -xml $output_directory/allNodeDisps.out -time -dof 1 2 3 4 5 6 disp]
# MODEL DETAILS RECORDER
print -JSON -file $output_directory/modelDetails.json
#
# #+# MODAL ANALYSIS BEFORE GRAVITY ------------------------------------------------
set nmodes 8; # Number of modes to analyze for modal analysis
# # Period before gravity
# if {$override != 2} {
# 	set wa [eigen $nmodes];
# 	set Periods 	[open $output_directory/PeriodsPreG.txt w];
# 	puts "\tFundamental-Period Before Gravity Analysis:"
# 	for {set iPd 1} {$iPd <= $nmodes} {incr iPd 1} {
# 		set wwa [lindex $wa $iPd-1];
# 		set Ta [expr 2*$pi/sqrt($wwa)];
# 		puts "\tPeriod$iPd= $Ta"
# 		puts $Periods "$Ta";
# 	}
# 	close $Periods;
# 	write_modes $output_directory/modesPreG.yaml $nmodes
# }
#
# #+# GRAVITY LOADING --------------------------------------------------------------
source getVertCosines.tcl
set vecxz301 "0 0 1";
pattern Plain 999 Linear {
	# Cap Beam Element Loads
	set CSDir "./Dimensions/CapCS/";  							# Directory containing cap beam cross section information
	set CSType "Cap"; 											# Cross section type is cap beam
	# Bents 2-11
	for {set ib 2} {$ib <= 11} {incr ib 1} {
		lassign [ReadMPR $CSDir $CSType $ib {}] A Iy Iz J; # Cross Section properties
		lassign [getVertCosines [expr 10000*$ib+10] $vecxz301] cosy cosz cosx curEle;
		eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
		lassign [getVertCosines [expr 10000*$ib+20] $vecxz301] cosy cosz cosx curEle;
		eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	}
	# Bent 12
	lassign [ReadMPR $CSDir $CSType 12 {}] A Iy Iz J; # Cross Section properties
	for {set eleCtr 0; set eleNum "120010 120020 120030 120040";} {$eleCtr<=3} {incr eleCtr} {
		lassign [getVertCosines [lindex $eleNum $eleCtr] $vecxz301] cosy cosz cosx curEle;
		eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	}
	# Bent 13 NE
	lassign [ReadMPR $CSDir $CSType 13 {}] A Iy Iz J; # Cross Section properties
	lassign [getVertCosines 130010 $vecxz301] cosy cosz cosx curEle;
	eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	lassign [getVertCosines 130020 $vecxz301] cosy cosz cosx curEle;
	eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	# Bent 14 NE
	lassign [ReadMPR $CSDir $CSType 14 {}] A Iy Iz J; # Cross Section properties
	lassign [getVertCosines 140010 $vecxz301] cosy cosz cosx curEle;
	eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	lassign [getVertCosines 140020 $vecxz301] cosy cosz cosx curEle;
	eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
	# Deck Element Loads
	set CSDir "./Dimensions/DeckCS/";  							# Directory containing deck cross section information
	set CSType "Deck"; 											# Cross section type is deck
	# NE Side
	for {set ib 1} {$ib <= 14} {incr ib 1} {
		for {set ibb 0} {$ibb <= 4} {incr ibb 1} {
			lassign [ReadMPR $CSDir $CSType $ib "-NodeNum $ibb"] A Iy Iz J; # Cross Section properties
			lassign [getVertCosines [expr 10*$ib+$ibb] $vecxz301] cosy cosz cosx curEle;
			eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
			#eleLoad -ele [expr 10*$ib+$ibb] -type -beamUniform 0 [expr -$wconc*$A]; # Uniformly distributed element load acting in vertical (local y) direction of element
		} 
	}
	for {set ib 12} {$ib <= 14} {incr ib 1} {
		for {set ibb 5} {$ibb <= 9} {incr ibb 1} {
			lassign [ReadMPR $CSDir $CSType $ib "-NodeNum $ibb"] A Iy Iz J; # Cross Section properties
			lassign [getVertCosines [expr 10*$ib+$ibb] $vecxz301] cosy cosz cosx curEle;
			eleLoad -ele $curEle -type -beamUniform [expr -$wconc*$A*$cosy] [expr -$wconc*$A*$cosz] [expr -$wconc*$A*$cosx];
			#eleLoad -ele [expr 10*$ib+$ibb] -type -beamUniform 0 [expr -$wconc*$A]; # Uniformly distributed element load acting in vertical (local y) direction of element
		} 
	}
};
#
# #+# STATIC (GRAVITY) ANALYSIS ----------------------------------------------------
wipeAnalysis
test NormDispIncr 1.0e-8 10 0;	
algorithm Newton;	
integrator LoadControl 0.1;
numberer Plain;
constraints Transformation;
system SparseGeneral;
analysis Static;
analyze 10;
write_displacements "$output_directory/dispsGrav.yaml"
# write_displacements "$output_directory/dispsGrav.yaml" Accel
#
# #+# MODAL ANALYSIS AFTER GRAVITY -------------------------------------------------
if {$override != 2} {
	# Period after gravity
	set wb [eigen -fullGenLapack $nmodes];
	set Periods 	[open $output_directory/PeriodsPostG.txt w];
	puts "\tFundamental-Period After Gravity Analysis:"
	for {set iPd 1} {$iPd <= $nmodes} {incr iPd 1} {
		set wwb [lindex $wb $iPd-1];
		set Tb [expr 2*$pi/sqrt($wwb)];
		puts "\tPeriod$iPd= $Tb"
		puts $Periods "$Tb";
	}
	close $Periods;
	write_modes $output_directory/modesPostG.yaml $nmodes
	for {set ctrRec 0} {$ctrRec<[llength $gravRec]} {incr ctrRec} {
		remove recorder [lindex $gravRec $ctrRec]
	}
}
#
# #+# TIMER 2 ----------------------------------------------------------------------
# Record time post gravity definition, static, and modal analysis
set tStatic [expr {([clock clicks -millisec]-$t1)/1000.}];
puts stderr "\tThe static and modal analysis time was $tStatic seconds"; # Print out the static and modal analysis time
puts $timers "\"tStatic\": $tStatic,";
set t2 [clock clicks -millisec];
#
# #+# TRANSIENT (DYNAMIC EARTHQUAKE GROUND MOTION ACCELERATION INPUT) ANALYSIS -----
if {$dynamic_on != 0} {
loadConst -time 0.0; # maintain constant gravity loads and reset time to zero
wipeAnalysis
# RAYLEIGH OR MODAL DAMPING
set nmodes [tcl::mathfunc::max {*}$damping_modes $nmodes]
# TODO: Consider removing fullGenLapack
set lambdaN [eigen -fullGenLapack $nmodes];
if {$damping_type == "rayleigh"} {
	set nEigenI [lindex $damping_modes 0]; 								# first rayleigh damping mode
	set nEigenJ [lindex $damping_modes 1]; 								# second rayleigh damping mode	
	set iDamp   [lindex $damping_ratios 0]; 							# first rayleigh damping ratio
	set jDamp   [lindex $damping_ratios 1]; 							# second rayleigh damping ratio
	set lambdaI [lindex $lambdaN [expr $nEigenI-1]];
	set lambdaJ [lindex $lambdaN [expr $nEigenJ-1]];
	set omegaI [expr $lambdaI**0.5];
	set omegaJ [expr $lambdaJ**0.5];
	set TI [expr 2.0*$pi/$omegaI];
	set TJ [expr 2.0*$pi/$omegaJ];
	set alpha0 [expr 2.0*($iDamp/$omegaI-$jDamp/$omegaJ)/(1/$omegaI**2-1/$omegaJ**2)];
	set alpha1 [expr 2.0*$iDamp/$omegaI-$alpha0/$omegaI**2];
	puts "\tRayleigh damping parameters:"
	puts "\tmodes: $nEigenI, $nEigenJ ; ratios: $iDamp, $jDamp"
	puts "\tTI = $TI; TJ = $TJ"
	puts "\tlambdaI = $lambdaI; lambdaJ = $lambdaJ"
	puts "\tomegaI = $omegaI; omegaJ = $omegaJ"
	puts "\talpha0 = $alpha0; alpha1 = $alpha1"
	rayleigh $alpha0 0.0 0.0 $alpha1;

} elseif {$damping_type == "modal"} {
	# needs a bit of edit. currently assuming that the ratios are applied in order at the first modes. but should be applied at the specified damping_modes modes.
	set nratios [llength $damping_ratios]
	puts "\tModal damping parameters:"
	puts "\tratios of $damping_ratios at the first $nratios modes"
	for {set i 1} {$i <= [expr $nmodes - $nratios]} {incr i} {
		lappend damping_ratios 0
	}
	modalDamping {*}$damping_ratios
}
#
# # #+# OUTPUT SYSTEM MATRICES (not yet implemented; matrices.tcl needs updating) ----------------------------------------------------------------------
# source matrices.tcl
# exit
# #
# #+# DYNAMIC ANALYSIS OPTIONS AND LOADING ------------------------------------------------------------
wipeAnalysis
set dtfact 1;
set SOS 1;
file mkdir $output_directory/model;
set Tol				1.0e-8;
set maxNumIter		100;
set printFlag		0;
set TestType		EnergyIncr;
set NewmarkGamma	0.50;
set NewmarkBeta		0.25;
constraints Transformation;
numberer RCM;
test $TestType $Tol $maxNumIter $printFlag;
if {$dynamic_integrator == "Newmark"} {
# set algorithmType   NewtonLineSearch;
set algorithmType   Newton;
# system SparseGeneral -piv;
system BandGeneral;
integrator Newmark $NewmarkGamma $NewmarkBeta;
} else {
	puts "\tspecify dynamic integrator"
};
# if {$dynamic_integrator == "OS"} {
# set algorithmType   Linear;
# system UmfPack;
# integrator AlphaOS 1.00;
# };
# if {$dynamic_integrator == "TR"} {
# set algorithmType   NewtonLineSearch;
# system SparseGeneral -piv;
# integrator TRBDF2;
# };
# if {$dynamic_integrator == "TR1"} {
# set algorithmType   NewtonLineSearch;
# system SparseGeneral -piv;
# integrator TRBDF3;
# };
# integrator HHT 0.7;
# integrator GeneralizedAlpha 1.0 0.8
algorithm $algorithmType;
analysis Transient;
if {$input_location != 0} {
	# Uniform Support Excitation
	lassign [py -m  CE58658.makePattern $record_zip --scale $dynamic_scale_factor --node $input_location] dt steps
} else {
	# # Multiple Support Excitation (not yet implemented; the below is a system ID impulse response experiment)
	# # IO1: same column and different directions
	# set exp changedof
	# set u1node 401
	# set u1dof 1
	# set u2node 401
	# set u2dof 2
	# set y1node 402
	# set y1dof 1
	# set y2node 402
	# set y2dof 2
	# set ui 1
	# set unode $u1node
	# set udof $u1dof
	# # set ui 2
	# # set unode $u2node
	# # set udof $u2dof
	# # # IO2: different column and same direction (Y)
	# # set exp changecol
	# # set u1node 401
	# # set u1dof 2
	# # set u2node 407
	# # set u2dof 2
	# # set y1node 402
	# # set y1dof 2
	# # set y2node 405
	# # set y2dof 2
	# # set ui 1
	# # set unode $u1node
	# # set udof $u1dof
	# # # set ui 2
	# # # set unode $u2node
	# # # set udof $u2dof
	# timeSeries Pulse 123456 0 1 10 -width 0.011 -factor 1
	# pattern MultipleSupport 123456 {
	# 	groundMotion 123456 Plain -accel 123456
	# 	imposedMotion $unode $udof 123456;
	# }
	# set dt 0.1
	# set steps 1000
	# algorithm $algorithmType;
	# analysis Transient;
	# set inFilelong 		"./Records/GM_long_global.txt";
	# set inFiletrans 	"./Records/GM_trans_global.txt";
	# set inFilevert 		"./Records/GM_vert_global.txt";
	# set iGMfile 		"$inFilelong $inFiletrans $inFilevert";
	# set iGMdirection 	"1 2 3";
	# set iGMfact 		"$SF $SF $SF";
	# set IDloadTag 		1;
	# set iloop 			"1 2 3";
	# set dt 				[expr 0.005*$sec];
	# set NumPts  		12600;
	# foreach Sloop $iloop GMdirection $iGMdirection GMfile $iGMfile GMfact $iGMfact {
	# 	incr IDloadTag;
	# 	set inFile	      "$GMfile.txt";
	# 	timeSeries Path $Sloop -dt $dt -filePath $GMfile -factor [expr $GMfact*$g];	
	# 	pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $Sloop;		# create Uniform excitation
	# 	if {$Sloop == 1} {
	# 		puts "\tHorizontal component 1 checked"
	# 	};
	# 	if {$Sloop == 2} {						
	# 		puts "\tHorizontal component 2 checked"
	# 	};
	# 	if {$Sloop == 3} {						
	# 		puts "\tVertical component checked"
	# 	};
	# };
	# set DtAnalysis 		[expr $dt/$dtfact];
	# set TmaxAnalysis 	[expr $dt*$NumPts];
	# set Nsteps 			[expr int($TmaxAnalysis/$DtAnalysis)];
	# puts "\tGround Motion: dt= $DtAnalysis, NumPts= [expr $NumPts*$dtfact], TmaxAnalysis= $TmaxAnalysis";
	# puts "\tRunning dynamic ground motion analysis..."
}
set DtAnalysis 	    $dt;
set TmaxAnalysis 	[expr $dt*$steps];
set Nsteps 			$steps;
if {$dynamic_truncated != 0} {
set Nsteps 			$dynamic_timesteps;
}
puts "\tGround Motion: dt= $DtAnalysis, NumPts= $Nsteps, TmaxAnalysis= $TmaxAnalysis";
#
# #+# DYNAMIC RECORDERS ------------------------------------------------------------
# TODO: Compute this dynamically
set column_nodes {203 205 303 305 403 405 503 505 603 605 703 705 803 805 903 905 1003 1005 1103 1105 1203 1207 1205 1303 1305 1315 1403 1405 1404 1415}
set column_elements {2010 2020 3010 3020 4010 4020 5010 5020 6010 6020 7010 7020 8010 8020 9010 9020 10010 10020 11010 11020 12010 12020 12030 13010 13020 13040 14010 14020 14030 14040}
# foreach node $column_nodes {
# 	catch {
# 		puts "NODE $node: [nodeCoord $node]"
# 	}
# }
set dynRec ""
## COLUMN 4010 DISPLACMENTS FOR PDCA DAMAGE STATES
# lappend dynRec [recorder Node -file $output_directory/Col4010TopDisp.txt -time -node 402 -dof 1 2 disp]
## COLUMN SECTION DEFORMATIONS AT TOP AND BOTTOM FOR STRAIN-BASED DAMAGE STATES
if {$column_linearity == "nonlinear"} {
	lappend dynRec [recorder Element -xml $output_directory/eleDef1.txt -ele {*}$column_elements section 1 deformation]
	lappend dynRec [recorder Element -xml $output_directory/eleDef4.txt -ele {*}$column_elements section 4 deformation]
}
# RESPONSE HISTORY RECORDERS
if {$input_location == 0} {
	# Impulse response experiment recorders
	file mkdir $output_directory/si_exp
	lappend dynRec [recorder Node -xml $output_directory/model/AA_all.txt -timeSeries 1 2 -dof 1 2 accel]
	lappend dynRec [recorder Node -xml $output_directory/model/RD_all.txt -dof 1 2 disp]
	lappend dynRec [recorder Node -file $output_directory/si_exp/u.txt -timeSeries 2 -node 401 -dof 2 accel]
	lappend dynRec [recorder Node -file $output_directory/si_exp/y1.txt -timeSeries 2 -node 402 -dof 2 accel]
	lappend dynRec [recorder Node -file $output_directory/si_exp/y2.txt -timeSeries 2 -node 405 -dof 2 accel]
} else {
	lappend dynRec [recorder Node -xml $output_directory/model/AA_all.txt -timeSeries 1 2 -dof 1 2 accel]
	# lappend dynRec [recorder Node -xml $output_directory/model/AA_all.txt -dof 1 2 accel]; # FOR DEBUGGING ONLY
	lappend dynRec [recorder Node -xml $output_directory/model/RD_all.txt -dof 1 2 disp]
}
# # Relative Displacement
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch02-03_X.txt -node 1031 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch02-03_Y.txt -node 1031 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch12-13_X.txt -node 1030 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch12-13_Y.txt -node 1030 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch06-07_X.txt -node 307 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch06-07_Y.txt -node 307 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch14-15_X.txt -node 304 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch14-15_Y.txt -node 304 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch17-18_X.txt -node 401 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch17-18_Y.txt -node 401 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch19-20_X.txt -node 402 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch19-20_Y.txt -node 402 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch24-25_X.txt -node 407 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch24-25_Y.txt -node 407 -dof 2 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch22-23_X.txt -node 405 -dof 1 disp]
# lappend dynRec [recorder Node -file $output_directory/model/RD_Ch22-23_Y.txt -node 405 -dof 2 disp]
# # Absolute Acceleration (with uniform excitation)
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch02-03_X.txt -timeSeries 1 -node 1031 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch02-03_Y.txt -timeSeries 2 -node 1031 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch12-13_X.txt -timeSeries 1 -node 1030 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch12-13_Y.txt -timeSeries 2 -node 1030 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch06-07_X.txt -timeSeries 1 -node 307 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch06-07_Y.txt -timeSeries 2 -node 307 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch14-15_X.txt -timeSeries 1 -node 304 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch14-15_Y.txt -timeSeries 2 -node 304 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch17-18_X.txt -timeSeries 1 -node 401 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch17-18_Y.txt -timeSeries 2 -node 401 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch19-20_X.txt -timeSeries 1 -node 402 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch19-20_Y.txt -timeSeries 2 -node 402 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch24-25_X.txt -timeSeries 1 -node 407 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch24-25_Y.txt -timeSeries 2 -node 407 -dof 2 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch22-23_X.txt -timeSeries 1 -node 405 -dof 1 accel]
# lappend dynRec [recorder Node -file $output_directory/model/AA_Ch22-23_Y.txt -timeSeries 2 -node 405 -dof 2 accel]
if {$input_location == 0} {
	# Impulse response experiment recorders
	file mkdir $output_directory/si_exp/$exp
	lappend dynRec [recorder Node -file $output_directory/si_exp/$exp/u1_u$ui.txt -node $u1node -dof $u1dof accel]
	lappend dynRec [recorder Node -file $output_directory/si_exp/$exp/u2_u$ui.txt -node $u2node -dof $u2dof accel]
	lappend dynRec [recorder Node -file $output_directory/si_exp/$exp/y1_u$ui.txt -node $y1node -dof $y1dof accel]
	lappend dynRec [recorder Node -file $output_directory/si_exp/$exp/y2_u$ui.txt -node $y2node -dof $y2dof accel]
}
# lappend dynRec [recorder Node -xml $output_directory/model/AA_X.txt -timeSeries 1 -node 1031 1030 307 304 401 402 407 405 -dof 1 accel]
# lappend dynRec [recorder Node -xml $output_directory/model/AA_Y.txt -timeSeries 2 -node 1031 1030 307 304 401 402 407 405 -dof 2 accel]
# TOP OF COLUMN ACCELERATION AND DRIFT RECORDERS
# lappend dynRec [recorder Node -xml  $output_directory/TopColAccel_X.txt -timeSeries 1 -node {*}$column_nodes -dof 1 accel]
# lappend dynRec [recorder Node -xml  $output_directory/TopColAccel_Y.txt -timeSeries 2 -node {*}$column_nodes -dof 2 accel]
# lappend dynRec [recorder Node -xml  $output_directory/TopColDrift_X.txt -node {*}$column_nodes -dof 1 disp]
# lappend dynRec [recorder Node -xml  $output_directory/TopColDrift_Y.txt -node {*}$column_nodes -dof 2 disp]
lappend dynRec [recorder Node -file $output_directory/TopColAccel_X_txt.txt -timeSeries 1 -node {*}$column_nodes -dof 1 accel]
lappend dynRec [recorder Node -file $output_directory/TopColAccel_Y_txt.txt -timeSeries 2 -node {*}$column_nodes -dof 2 accel]
lappend dynRec [recorder Node -file $output_directory/TopColDrift_X_txt.txt -node {*}$column_nodes -dof 1 disp]
lappend dynRec [recorder Node -file $output_directory/TopColDrift_Y_txt.txt -node {*}$column_nodes -dof 2 disp]
# ALL NODES BINARY RECORDER (large file! be careful)
# lappend dynRec [recorder Node -binary $output_directory/allNodes.bin		-precision 5	-time	-dof 1 2 3 4 5 6 disp]
#
# #+# TIMER 3 ----------------------------------------------------------------------
# Record time for Rayleigh damping and dynamic analysis definition;
set tDynDef [expr {([clock clicks -millisec]-$t2)/1000.}];
puts stderr "\tRayleigh or modal damping and dynamic analysis definition time was $tDynDef seconds";
puts $timers "\"tDynDef\": $tDynDef,";
#
# #+# SOLVE DYNAMIC ANALYSIS -------------------------------------------------------
# Start timer for the dynamic analysis
puts "\tRunning dynamic ground motion analysis..."
set t3 [clock clicks -millisec];
catch {progress create $Nsteps} _
for {set ik 1} {$ik <= $Nsteps} {incr ik 1} {
	catch {progress update} _ 
      set ok	  [analyze 1 $DtAnalysis];
	# Convergence
    if {$ok != 0} {
		puts "\tTrying Bisection ...";
		algorithm NewtonLineSearch -type Bisection;
		set ok [analyze 1 $DtAnalysis]
	};
	if {$ok != 0} {
		puts "\tTrying Secant ...";
		algorithm NewtonLineSearch -type Secant;
		set ok [analyze 1 $DtAnalysis]
	};
	if {$ok != 0} {
		puts "\tTrying RegulaFalsi ...";
		algorithm NewtonLineSearch -type RegulaFalsi;
		set ok [analyze 1 $DtAnalysis]
	};
	if {$ok != 0} {
		puts "\tTrying KrylovNewton ...";
		algorithm KrylovNewton;
		set ok [analyze 1 $DtAnalysis];
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "\tTrying Newton...";
		algorithm Newton;
		set ok [analyze 1 $DtAnalysis];
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "\tTrying Broyden ...";
		algorithm Broyden;
		set ok [analyze 1 $DtAnalysis];
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "\tTrying BFGS ...";
		algorithm BFGS;
		set ok [analyze 1 $DtAnalysis];
		algorithm $algorithmType;
	};
	if {$SOS == 1} {
	    if {$ok != 0} {
		    puts "\tTrying OS ...";
			integrator AlphaOS 1.00;
		    algorithm Linear;
		    set ok [analyze 1 $DtAnalysis];
			integrator Newmark $NewmarkGamma $NewmarkBeta;
		    algorithm $algorithmType;
	    };
	};
	if {$ok != 0} {
		set Nstepsmax [expr $ik-1];
		break;
	};
};
if {[expr $ik-1] == $Nsteps} {
	set AnalysisA [expr 1] 
	puts "\tAnalysis complete"
} else {
	set AnalysisA [expr 0]
};
#
# #+# MODAL ANALYSIS AFTER DYNAMIC ANALYSIS -------------------------------------------------
# # Period after dynamic analysis
# set wb [eigen -fullGenLapack $nmodes];
# set Periods 	[open $output_directory/PeriodsPostD.txt w];
# puts "\tFundamental-Period After Dynamic Analysis:"
# for {set iPd 1} {$iPd <= $nmodes} {incr iPd 1} {
# 	set wwb [lindex $wb $iPd-1];
# 	set Tb [expr 2*$pi/sqrt($wwb)];
# 	puts "\tPeriod$iPd= $Tb"
# 	puts $Periods "$Tb";
# }
# close $Periods;
# write_modes $output_directory/modesPostD.yaml $nmodes
#
# #+# TIMER 4 ----------------------------------------------------------------------
# Record time post dynamic analysis
set tDynamic [expr {([clock clicks -millisec]-$t3)/1000.}];
puts stderr "\tThe dynamic analysis time was $tDynamic seconds"; # Print out the dynamic analysis time
puts $timers "\"tDynamic\": $tDynamic,";
};
# Record total time for all model and analysis definition and execution
set tTotal [expr {([clock clicks -millisec]-$t0)/1000.}];
puts stderr "\tThe total time for all model and analysis definition and execution was $tTotal seconds"; # Print out the total time
puts $timers "\"tTotal\": $tTotal}";
close $timers;
#
# #+# END ----------------------------------------------------------------------
#
wipe