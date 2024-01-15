	#########################################################
	# Displacement history file for one cycle of displacement
	#########################################################
	
	# Written: Vesna Terzic (vesna@berkeley.edu)
	# Created: 12/2011

proc singlecycle {umax n running} {

	# note, m should always equal the rate, by definition
	set rate 0.01
	set N [expr 4*$n-3]
	set dt [expr 4.0*$umax/$rate/($N-1)]
	set m [expr $umax/($n-1)/$dt]

	set filet [open out/time.txt "a"]
	set fileu [open out/displacement.txt "a"]

	for { set i 1 } { $i <= $N } { incr i } {
		set tval($i) [expr ($i-1)*$dt]
		if { $i >= [expr 3*$n-2] } {
			set uval($i) [expr $m*$tval($i)-4.0*$umax]
		} elseif { $i >= [expr $n+1]} {
			set uval($i) [expr -$m*$tval($i)+2.0*$umax]
		} else {
			set uval($i) [expr $m*$tval($i)]
		}
		set tout($i) $tval($i)
		set tval($i) [expr $tval($i) + $running]
		puts $fileu $uval($i)
		puts $filet $tval($i)
	}

	close $fileu
	close $filet

	return $tout($N)

}
