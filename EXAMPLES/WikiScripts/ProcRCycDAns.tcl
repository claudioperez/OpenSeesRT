#####################################################################################################
# #
# procRCycDAns.tcl #
# procedure for reverse cyclic displacement control analysis given the peak pts. #
# analysis type used : STATIC #
# Written : N.Mitra #
#####################################################################################################
proc procRCycDAns { incre nodeTag dofTag peakpts} {
set x [lindex $peakpts 0]
set fir [expr $x/$incre]
integrator DisplacementControl $nodeTag $dofTag 0.0 1 $fir $fir
# create the analysis object
analysis Static
# perform the analysis
analyze $incre
integrator DisplacementControl $nodeTag $dofTag 0.0 1 [expr -$fir] [expr -$fir]
analyze [expr 2*$incre]
integrator DisplacementControl $nodeTag $dofTag 0.0 1 $fir $fir
analyze $incre
for {set j 1} {$j < [llength $peakpts]} {incr j 1} {
set tx [lindex $peakpts $j]
set tinc [expr $tx/$fir]
set rt [expr int($tinc)]
integrator DisplacementControl $nodeTag $dofTag 0.0 1 $fir $fir
analyze $rt
integrator DisplacementControl $nodeTag $dofTag 0.0 1 [expr -$fir] [expr -$fir]
analyze [expr 2*$rt]
integrator DisplacementControl $nodeTag $dofTag 0.0 1 $fir $fir
analyze $rt
}
################################ end procRCycDAns.tcl #######################################
}