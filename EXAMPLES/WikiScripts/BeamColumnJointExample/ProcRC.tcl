########################################################################################################

#

# procRC.tcl

## procedure for setting up a reversed cycle loading scheme. The input are mainly the

## peak points for the loading.

## The procedure primarily uses Displacement control for loading, if it fails uses ArcLength control

## created : NM (nmitra@u.washington.edu) dated: Sep 2002

########################################################################################################

proc procRC { incre nodeTag dofTag peakpts } {

set displayTag 0;

set numTimes 150;

set x [lindex $peakpts 0];

set dU [expr $x/$incre];

#set dU0 [expr $dU/1000];

set dU0 [expr $dU/10000];

integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU $dU

analysis Static

analyze $incre

integrator DisplacementControl $nodeTag $dofTag 0.0 1 [expr -$dU] [expr -$dU]

analyze [expr 2*$incre]

integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU $dU

analyze $incre

## end the first peak pt start for others

for {set j 1} {$j < [llength $peakpts]} {incr j 1} {

set y [lindex $peakpts $j]

set dSt [expr $y/$dU]

set dS [expr int($dSt)]

test NormDispIncr 1e-8 $numTimes $displayTag

algorithm Newton

############# start loading cycle ##################

set t 0;

while {$t != $dS} {

integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU $dU

set ok [analyze 1]

incr t 1;

if {$ok != 0} {

# if {$t == $dS} {break};

puts "Displacement control failed ..... trying Arc-Length control"

set currentDisp [nodeDisp $nodeTag $dofTag]

puts "Current Displacement is $currentDisp"

# algorithm Linear

test NormDispIncr 1e-6 $numTimes $displayTag

#algorithm ModifiedNewton

# integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU0 $dU0

# integrator DisplacementControl $nodeTag $dofTag 0.0 10 $dU0 $dU0

integrator ArcLength [expr $dU0] 1.0

# set ok [analyze 1]

analyze 1

}

# puts "that worked ..... back to regular Newton "

test NormDispIncr 1e-8 $numTimes $displayTag

# algorithm Newton

}

################## end of loading cycle, start unloading cycle ########

set t 0;

while {$t != [expr 2*$dS]} {

integrator DisplacementControl $nodeTag $dofTag 0.0 1 [expr -$dU] [expr -$dU]

set ok [analyze 1]

incr t 1;

if {$ok != 0} {

# if {$t == [expr 2*$dS]} {break};

puts "Displacement control failed ..... trying Arc-Length control"

set currentDisp [nodeDisp $nodeTag $dofTag]

puts "Current Displacement is $currentDisp"

# algorithm Linear

test NormDispIncr 1e-6 $numTimes $displayTag

#algorithm ModifiedNewton

# integrator DisplacementControl $nodeTag $dofTag 0.0 1 [expr -$dU0] [expr -$dU0]

# integrator DisplacementControl $nodeTag $dofTag 0.0 10 [expr -$dU0] [expr -$dU0]

integrator ArcLength [expr $dU0] 1.0

# set ok [analyze 1]

analyze 1

}

# puts "that worked .... back to regular Newton "

test NormDispIncr 1e-8 $numTimes $displayTag

# algorithm Newton

}

############# end of unloading cycle, start reloading cycle ###########

set t 0;

while {$t != $dS} {

integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU $dU

set ok [analyze 1]

incr t 1;

if {$ok != 0} {

# if {$t == $dS} {break};

puts "Displacement control failed ..... trying Arc-Length control"

set currentDisp [nodeDisp $nodeTag $dofTag]

puts "Current Displacement is $currentDisp"

# algorithm Linear

test NormDispIncr 1e-6 $numTimes $displayTag

#algorithm ModifiedNewton

# integrator DisplacementControl $nodeTag $dofTag 0.0 1 $dU0 $dU0

# integrator DisplacementControl $nodeTag $dofTag 0.0 10 $dU0 $dU0

integrator ArcLength [expr $dU0] 1.0

# set ok [analyze 1]

analyze 1

}

# puts "that worked .... back to regular Newton "

test NormDispIncr 1e-8 $numTimes $displayTag

# algorithm Newton

}

######## reloading cycle completed #############################

if {$ok == 0} {

puts "analysis succesful at $y mm displacement";

} else {

puts "analysis could not proceed fine beyond $y mm displacement";

}

}

}