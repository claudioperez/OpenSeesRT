set length 500.
set height 100

set area  10.
set young 100000.
set mat 1

set nx 10000
set ny 2

set Px 5
set Py 5

set k 0
set tag 0
set truss "$area $mat"

model basic -ndm 2 -ndf 2
uniaxialMaterial Elastic $mat $young

node [incr k] 0.0 0.0
node [incr k] 0.0 $height
fix 1 1 1
fix 2 1 1

element Truss [incr tag] [expr $k-1] [expr $k] {*}$truss ; # vertical

for {set i 1} {$i < $nx} {incr i} {
  node [incr k] [expr $i*$length/$nx] 0.0
  node [incr k] [expr $i*$length/$nx] $height
  
  element Truss [incr tag] [expr $k-3] [expr $k] {*}$truss   ; # diagonal
  element Truss [incr tag] [expr $k-1] [expr $k] {*}$truss   ; # vertical
  element Truss [incr tag] [expr $k-2] [expr $k] {*}$truss   ; # horizontal
  element Truss [incr tag] [expr $k-1] [expr $k-3] {*}$truss ; # horizontal
}

print -JSON test.json


# Define loads
# ------------
timeSeries Linear 1
pattern Plain 1 1 {
   # Create the nodal load - command: load nodeID xForce yForce
   load $k $Px -$Py
}


system petsc
numberer RCM

integrator LoadControl 1.0
analysis Static

start
analyze 1
stop


# foreach node [getNodeTags] {
#   puts "$node: \[[join [nodeDisp $node] {, } ] , 0.0, 0.0, 0.0, 0.0 \]"
# }


