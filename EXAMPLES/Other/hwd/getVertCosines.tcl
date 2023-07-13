proc getVertCosines {eleTag  vecxz} {
	set nodeTags [eleNodes $eleTag]
	set n1Crds	[nodeCoord [lindex $nodeTags 0]]
	set n2Crds	[nodeCoord [lindex $nodeTags 1]]
	set eleLocX(0)	[expr [lindex $n2Crds 0]-[lindex $n1Crds 0]]
	set eleLocX(1)	[expr [lindex $n2Crds 1]-[lindex $n1Crds 1]]
	set eleLocX(2)	[expr [lindex $n2Crds 2]-[lindex $n1Crds 2]]
	set eleLocXMag [expr sqrt($eleLocX(0)**2+$eleLocX(1)**2+$eleLocX(2)**2)]
	set eleLocXNorm(0)	[expr $eleLocX(0)/$eleLocXMag]
	set eleLocXNorm(1)	[expr $eleLocX(1)/$eleLocXMag]
	set eleLocXNorm(2)	[expr $eleLocX(2)/$eleLocXMag]
	
	set vecxzMag [expr sqrt([lindex $vecxz 0]**2+[lindex $vecxz 1]**2+[lindex $vecxz 2]**2)]

	array set vecxzNorm "
		0	[expr [lindex $vecxz 0]/$vecxzMag]
		1	[expr [lindex $vecxz 1]/$vecxzMag]
		2	[expr [lindex $vecxz 2]/$vecxzMag]
	"
	array set eleLocYNorm "
		0	[expr       $vecxzNorm(1)*$eleLocXNorm(2)-$vecxzNorm(2)*$eleLocXNorm(1)  ]
		1	[expr -1 * ($vecxzNorm(0)*$eleLocXNorm(2)-$vecxzNorm(2)*$eleLocXNorm(0)) ]
		2	[expr       $vecxzNorm(0)*$eleLocXNorm(1)-$vecxzNorm(1)*$eleLocXNorm(0)  ]
	"
	array set eleLocZNorm "
		0	[expr -1 * ($eleLocYNorm(1)*$eleLocXNorm(2)-$eleLocYNorm(2)*$eleLocXNorm(1)) ]
		1	[expr  1 * ($eleLocYNorm(0)*$eleLocXNorm(2)-$eleLocYNorm(2)*$eleLocXNorm(0)) ]
		2	[expr -1 * ($eleLocYNorm(0)*$eleLocXNorm(1)-$eleLocYNorm(1)*$eleLocXNorm(0)) ]
	"
#	puts "$eleLocXNorm(0) $eleLocXNorm(1) $eleLocXNorm(2)"
#	puts "$eleLocYNorm(0) $eleLocYNorm(1) $eleLocYNorm(2)"
#	puts "$eleLocZNorm(0) $eleLocZNorm(1) $eleLocZNorm(2)"
	return "$eleLocYNorm(2) $eleLocZNorm(2) $eleLocXNorm(2) $eleTag $eleLocXMag"
}
