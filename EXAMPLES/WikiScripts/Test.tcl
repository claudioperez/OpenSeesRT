## Units kips-inches
##----------------------------------------------------------
## Source in the model
source SAWSZeroLength.tcl
## Apply the nodal Load
pattern Plain 1 Linear { load 2 1.0 0.0 0.0 }
## Recorder
recorder Node -file out/LoadDisp.dat -time -node 2 -dof 1 disp
## Static Analysis parameters
test EnergyIncr 1.0e-8   300    0
algorithm KrylovNewton
system UmfPack
numberer RCM
constraints Plain
analysis Static
set peaks [list 0.1] 
for {set cntr 1} { $cntr < 20 } {incr cntr} {
    for {set cntr2 1} { $cntr2 <= 6 } {incr cntr2} {
	set peaks [linsert $peaks $cntr [expr $cntr*0.1 ]  ]
    }
}
for {set i 1 } { $i <= [llength $peaks] } {incr i } {
    set dU [expr -1.0*pow((-1.0),$i)*[lindex $peaks [expr $i-1] ]/50.0]
    integrator DisplacementControl 2 1 $dU 1 $dU $dU
    analyze 50
}
