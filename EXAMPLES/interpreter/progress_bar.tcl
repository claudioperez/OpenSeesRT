
set delay 100 ; # 0.1 seconds

# progress create
# for {set i 0} {$i < 100} {incr i} {after $delay; progress update}
# puts ""

progress create 200
for {set i 0} {$i <  50} {incr i} {after $delay; progress update ""}
for {set i 0} {$i <  50} {incr i} {after $delay; progress update "Part I"}
for {set i 0} {$i <  50} {incr i} {after $delay; progress update "Part II"}
for {set i 0} {$i <  50} {incr i} {after $delay; progress update ""}

puts ""
