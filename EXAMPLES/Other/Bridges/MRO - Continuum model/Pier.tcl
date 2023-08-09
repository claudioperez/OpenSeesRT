
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates the bridge pier per site condition at Meloland Road Overpass, in El Centro, CA.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================

################################
###### Develop Pier Nodes ######
################################

# Creates input file for GiD software
puts $meshFile "MESH beamcolumn dimension 3 ElemType linear Nnode 2"
puts $meshFile "Coordinates"


set nodeNum22  $nodeNum176
set nodeNum24   $nodeNum22 

for { set i 1 } { $i <= [ expr $Hpier/$psize+$Df/$zSize3+1] } { incr i 1 } {

    if { $i < [expr $Df/$zSize3+1] } { 
    node   $nodeNum24 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile + $WInfpile+$Lzone11/2.0] [ expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11/2.0]  [expr $HZ-$Df+($i-1)*$zSize3]
    puts  $meshFile "$nodeNum24 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile + $WInfpile+$Lzone11/2.0] [ expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11/2.0]  [expr $HZ-$Df+($i-1)*$zSize3]"
    fix   $nodeNum24 0 0 0 0 0 1
    equalDOF  [expr $pierbase+($i-1)]  $nodeNum24 1 2 3
    set   nodeNum24    [expr $nodeNum24+1]
    } else {
    node  $nodeNum24 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile + $WInfpile+$Lzone11/2.0] [ expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11/2.0]  [expr $HZ+($i-($Df/$zSize3+1))*$psize]
    puts  $meshFile "$nodeNum24 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile + $WInfpile+$Lzone11/2.0] [ expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11/2.0]  [expr $HZ+($i-($Df/$zSize3+1))*$psize]"
    fix   $nodeNum24  0 0 0 0 0 1
    set   nodeNum24    [expr $nodeNum24+1]
    }
 }

puts $meshFile "end coordinates"


###################################
###### Develop Pier Elements ######
###################################

set nepier $nodeNum22
set eleNum [expr $eleNum+1]
set rec_for_pier $eleNum
puts $meshFile "Elements"

for { set i 1 } { $i <= [ expr $Hpier/$psize+$Df/$zSize3] } { incr i 1 } {

  element dispBeamColumn $eleNum $nepier [expr $nepier+1]   $numIntgrPts 2  101  -mass $massDens_pier 
  puts $meshFile "$eleNum $nepier [expr $nepier+1] 10"

  set nepier [expr $nepier+1]
  set eleNum [expr $eleNum+1]
    }

puts $meshFile "end elements"     
  



