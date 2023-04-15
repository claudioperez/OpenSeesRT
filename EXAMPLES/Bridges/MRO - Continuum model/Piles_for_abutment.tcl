
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates pile groups under left and right abutment walls per site condition at Meloland Road Overpass, in El Centro, CA.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================

#############################################
###### Develop Pile Nodes at Abutments ######
#############################################

# Creates input file for GiD software
puts $meshFile "MESH beamcolumn dimension 3 ElemType linear Nnode 2"
puts $meshFile "Coordinates"

set pile_elenum_abut [expr int($Df_abut/$zSize4+$Czone3/$zSize3+$Czone2/$zSize2+($LP_abut-$Czone3-$Czone2)/$zSize1)]

# See my Ph.D. dissertation for the pile labeling

# Left Abutment Piles	
########################################################
#################### PILE A ############################
########################################################

set   nodeNum149    [expr $nodeNum148+1]

set m 1
set n $p1_a
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a }
  if { $i == 3 } { set n $p3_a }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+($j-1)*$xSize4 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum149 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum149 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a $nodeNum149       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum149 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum149 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum149 1 2 3 
            } 
         set nodeNum149 [expr $nodeNum149+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+(2-1)*$xSize4 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum149  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum149 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+1)*$numZele4+$m)]   $nodeNum149    1  2  3 
         set m [expr $m+1]
         set nodeNum149 [expr $nodeNum149+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum150    $nodeNum149

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum150 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+$BP_emb/2.0]  $zdim
      equalDOF  $nodeNum150  [expr $np_a+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum150 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+$BP_emb/2.0]  $zdim"
      set  nodeNum150  [expr $nodeNum150+1]
    } 


for { set i $np_a } { $i <= [expr int($np_a+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE A1 ############################
########################################################

set   nodeNum151    $nodeNum150

set m 1
set n $p1_a1
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a1 }
  if { $i == 3 } { set n $p3_a1 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum151 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum151 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a1 $nodeNum151       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum151 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum151 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum151 1 2 3 
            } 
         set nodeNum151 [expr $nodeNum151+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum151  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum151 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+4+1)*$numZele4+$m)]   $nodeNum151    1  2  3 
         set m [expr $m+1]
         set nodeNum151 [expr $nodeNum151+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum152    $nodeNum151

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum152 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum152  [expr $np_a1+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum152 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum152  [expr $nodeNum152+1]
    } 


for { set i $np_a1 } { $i <= [expr int($np_a1+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 
	
########################################################
#################### PILE A2 ############################
########################################################

set   nodeNum153    $nodeNum152

set m 1
set n $p1_a2
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a2 }
  if { $i == 3 } { set n $p3_a2 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum153 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum153 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a2 $nodeNum153       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum153 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum153 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum153 1 2 3 
            } 
         set nodeNum153 [expr $nodeNum153+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum153  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum153 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+8+1)*$numZele4+$m)]   $nodeNum153    1  2  3 
         set m [expr $m+1]
         set nodeNum153 [expr $nodeNum153+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum154    $nodeNum153

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum154 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum154  [expr $np_a2+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum154 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum154  [expr $nodeNum154+1]
    } 


for { set i $np_a2 } { $i <= [expr int($np_a2+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE A3 ############################
########################################################

set   nodeNum155    $nodeNum154

set m 1
set n $p1_a3
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a3 }
  if { $i == 3 } { set n $p3_a3 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]


         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum155 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum155 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a3 $nodeNum155       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum155 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum155 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum155 1 2 3 
            } 
         set nodeNum155 [expr $nodeNum155+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum155  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum155 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+11+1)*$numZele4+$m)]   $nodeNum155    1  2  3 
         set m [expr $m+1]
         set nodeNum155 [expr $nodeNum155+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum156    $nodeNum155

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum156 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum156  [expr $np_a3+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum156 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum156  [expr $nodeNum156+1]
    } 


for { set i $np_a3 } { $i <= [expr int($np_a3+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE A4 ############################
########################################################

set   nodeNum157    $nodeNum156

set m 1
set n $p1_a4
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a4 }
  if { $i == 3 } { set n $p3_a4 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum157 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum157 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a4 $nodeNum157      
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum157 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum157 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum157 1 2 3 
            } 
         set nodeNum157 [expr $nodeNum157+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum157  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum157 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+14+1)*$numZele4+$m)]   $nodeNum157    1  2  3 
         set m [expr $m+1]
         set nodeNum157 [expr $nodeNum157+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum158    $nodeNum157

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum158 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum158  [expr $np_a4+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum158 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum158  [expr $nodeNum158+1]
    } 


for { set i $np_a4 } { $i <= [expr int($np_a4+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 
	
########################################################
#################### PILE A5 ############################
########################################################

set   nodeNum159    $nodeNum158

set m 1
set n $p1_a5
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_a5 }
  if { $i == 3 } { set n $p3_a5 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum159 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum159 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_a5 $nodeNum159       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum159 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum159 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum159 1 2 3 
            } 
         set nodeNum159 [expr $nodeNum159+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum159  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum159 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+18+1)*$numZele4+$m)]   $nodeNum159    1  2  3 
         set m [expr $m+1]
         set nodeNum159 [expr $nodeNum159+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum160    $nodeNum159

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum160 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ]   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      equalDOF  $nodeNum160  [expr $np_a5+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum160 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim"
      set  nodeNum160  [expr $nodeNum160+1]
    } 


for { set i $np_a5 } { $i <= [expr int($np_a5+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 	

########################################################
#################### PILE B ############################
########################################################

set   nodeNum161    $nodeNum160
set m 1
set n $p1_b
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_b }
  if { $i == 3 } { set n $p3_b }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
 
      set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-1)*$xSize4 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+($j-1)*$xSize4 ]

      if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum161 $xdim $ydim $zdim
      puts $meshFile "$nodeNum161 $xdim $ydim $zdim"
      
      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_b $nodeNum161       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum161 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum161 1 2 3 
         }  
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum161 1 2 3 
            } 
       set nodeNum161 [expr $nodeNum161+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+(2-1)*$xSize4 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum161  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum161 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum1-1)+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($numYele3+7+$numXele11+5+1)*$numZele4+$m)]   $nodeNum161    1  2  3 
         set m [expr $m+1]
         set nodeNum161 [expr $nodeNum161+1]
       }}


    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum162    $nodeNum161

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum162 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+ $BP_emb/2.0 ]  $zdim
      equalDOF  $nodeNum162  [expr $np_b+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum162 [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0 ] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+ $BP_emb/2.0 ]  $zdim"
      set  nodeNum162  [expr $nodeNum162+1]
    } 


for { set i $np_b } { $i <= [expr int($np_b+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 


# Right Abutment Piles	
########################################################
#################### PILE G ############################
########################################################
set   nodeNum163   $nodeNum162

set m 1
set n $p1_g
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g }
  if { $i == 3 } { set n $p3_g }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {
 
      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

      set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+($j-1)*$xSize4 ]

      if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum163 $xdim $ydim $zdim
      puts $meshFile "$nodeNum163 $xdim $ydim $zdim"
      
      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_g $nodeNum163       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum163 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum163 1 2 3 
         }  
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum163 1 2 3 
            } 
       set nodeNum163 [expr $nodeNum163+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+(2-1)*$xSize4 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum163  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum163 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+1)*$numZele4+$m)]   $nodeNum163    1  2  3 
         set m [expr $m+1]
         set nodeNum163 [expr $nodeNum163+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum164    $nodeNum163

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum164 $L3 [expr $Bzone1+$Bzone2+$Bzone3+$BP_emb/2.0]  $zdim
      equalDOF  $nodeNum164  [expr $np_g+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum164 $L3 [expr $Bzone1+$Bzone2+$Bzone3+$BP_emb/2.0]  $zdim"
      set  nodeNum164  [expr $nodeNum164+1]
    } 


for { set i $np_g } { $i <= [expr int($np_g+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

	
########################################################
#################### PILE G1 ############################
########################################################

set   nodeNum165    $nodeNum164

set m 1
set n $p1_g1
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g1 }
  if { $i == 3 } { set n $p3_g1 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]
		 

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum165 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum165 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_g1 $nodeNum165       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum165 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum165 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum165 1 2 3 
            } 
         set nodeNum165 [expr $nodeNum165+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
        
		 set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum165  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum165 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+4+1)*$numZele4+$m)]   $nodeNum165    1  2  3 
         set m [expr $m+1]
         set nodeNum165 [expr $nodeNum165+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum166    $nodeNum165

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum166 $L3   [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum166  [expr $np_g1+$r-1]  1 2 3 

      puts $meshFile "$nodeNum166 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ] $zdim"
      set  nodeNum166  [expr $nodeNum166+1]
    } 


for { set i $np_g1 } { $i <= [expr int($np_g1+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 
	
########################################################
#################### PILE G2 ############################
########################################################

set   nodeNum167    $nodeNum166

set m 1
set n $p1_g2
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g2 }
  if { $i == 3 } { set n $p3_g2 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum167 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum167 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_g2 $nodeNum167       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum167 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum167 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum167 1 2 3 
            } 
         set nodeNum167 [expr $nodeNum167+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum167  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum167 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+8+1)*$numZele4+$m)]    $nodeNum167    1  2  3 
         set m [expr $m+1]
         set nodeNum167 [expr $nodeNum167+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum168    $nodeNum167

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum168 $L3  [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum168  [expr $np_g2+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum168 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum168  [expr $nodeNum168+1]
    } 


for { set i $np_g2 } { $i <= [expr int($np_g2+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE G3 ############################
########################################################

set   nodeNum169    $nodeNum168

set m 1
set n $p1_g3
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g3 }
  if { $i == 3 } { set n $p3_g3 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]


         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum169 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum169 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_g3 $nodeNum169       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum169 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum169 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum169 1 2 3 
            } 
         set nodeNum169 [expr $nodeNum169+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum169  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum169 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+11+1)*$numZele4+$m)]    $nodeNum169    1  2  3 
         set m [expr $m+1]
         set nodeNum169 [expr $nodeNum169+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum170    $nodeNum169

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum170 $L3  [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum170  [expr $np_g3+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum170 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum170  [expr $nodeNum170+1]
    } 


for { set i $np_g3 } { $i <= [expr int($np_g3+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE G4 ############################
########################################################

set   nodeNum171    $nodeNum170

set m 1
set n $p1_g4
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g4 }
  if { $i == 3 } { set n $p3_g4 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum171 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum171 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_g4 $nodeNum171      
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum171 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum171 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum171 1 2 3 
            } 
         set nodeNum171 [expr $nodeNum171+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum171  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum171 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+14+1)*$numZele4+$m)]   $nodeNum171    1  2  3 
         set m [expr $m+1]
         set nodeNum171 [expr $nodeNum171+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum172    $nodeNum171

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum172 $L3  [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      equalDOF  $nodeNum172  [expr $np_g4+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum172 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum172  [expr $nodeNum172+1]
    } 


for { set i $np_g4 } { $i <= [expr int($np_g4+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 
	
########################################################
#################### PILE G5 ############################
########################################################

set   nodeNum173    $nodeNum172

set m 1
set n $p1_g5
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_g5 }
  if { $i == 3 } { set n $p3_g5 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {

      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

         set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

         if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
            set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
         } else {
         if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
            set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
         } else {
            set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
         }}
      
         node $nodeNum173 $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum173 $xdim $ydim $zdim"
          
         if { $k == 1 && $i == 2 && $j == 2 } {
            set np_g5 $nodeNum173       
            }

         #######  Slaving Pile Nodes to Soil Nodes
         if { $i != 2 || $j != 2 } {
             equalDOF  [expr $n+$k-1] $nodeNum173 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == 1 } {
             equalDOF  [expr $n+$k-1] $nodeNum173 1 2 3 
            } 
         if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
             equalDOF  [expr $n+$k-1] $nodeNum173 1 2 3 
            } 
         set nodeNum173 [expr $nodeNum173+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
         set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+(2-1)*$xSize10 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum173  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum173 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+18+1)*$numZele4+$m)]    $nodeNum173    1  2  3 
         set m [expr $m+1]
         set nodeNum173 [expr $nodeNum173+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}


set nodeNum174    $nodeNum173

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum174 $L3  [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      equalDOF  $nodeNum174  [expr $np_g5+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum174 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim"
      set  nodeNum174  [expr $nodeNum174+1]
    } 


for { set i $np_g5 } { $i <= [expr int($np_g5+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 	
		
	
########################################################
#################### PILE H ############################
########################################################

set   nodeNum175    $nodeNum174

set m 1
set n $p1_h
for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_h }
  if { $i == 3 } { set n $p3_h }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_abut+1]} {incr k 1} {
 
      if { $k <= [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {

      set xdim [expr $L3-$BP_emb/2.0+($i-1)*$xSize4 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+($j-1)*$xSize4 ]

      if { $k <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $k > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum175 $xdim $ydim $zdim
      puts $meshFile "$nodeNum175 $xdim $ydim $zdim"
      
      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_h $nodeNum175       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum175 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum175 1 2 3 
         }  
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum175 1 2 3 
         } 
       set nodeNum175 [expr $nodeNum175+1]
       
       } else {
       if { $i == 2 && $j == 2 && $k > [expr $pile_elenum_abut-$Df_abut/$zSize4+1] } {
      set xdim [expr $L3-$BP_emb/2.0+(2-1)*$xSize4 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+(2-1)*$xSize4 ]
         set zdim [expr $HZ+($k-($pile_elenum_abut-$Df_abut/$zSize4+1))*$zSize4]

         node $nodeNum175  $xdim   $ydim   $zdim
         puts $meshFile "$nodeNum175 $xdim $ydim $zdim"

         equalDOF  [expr int(($nodeNum4-1)+($numYele3+7+$numXele11+5+1)*$numZele4+$m)]   $nodeNum175   1  2  3 
         set m [expr $m+1]
         set nodeNum175 [expr $nodeNum175+1]
       }}

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum176    $nodeNum175

for { set r 1 } { $r <= [expr ($pile_elenum_abut-$Df_abut/$zSize4+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   $LP_abut]
      } else {
        if { $r > [expr ($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_abut-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum176 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+ $BP_emb/2.0 ]  $zdim
      equalDOF  $nodeNum176  [expr $np_h+$r-1]  1 2 3 
 
      puts $meshFile "$nodeNum176 $L3 [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$Lzone11+2.0*$WInfemb+$BP_pile+$BP_pile/2.0+ $BP_emb/2.0 ]  $zdim"
      set  nodeNum176  [expr $nodeNum176+1]
    } 


for { set i $np_h } { $i <= [expr int($np_h+$pile_elenum_abut)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

puts $meshFile "end coordinates"

######################################################
###### Removing Soil Elements for Pile A #############
######################################################

set  numXELE_A1    [expr $numXele1+$numXele2+$numXele3+1]
set  numXELE_A2    $num_pile_Y1
for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile A1 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+3+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile A2 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+7+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile A3 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+10+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile A4 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+13+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile A5 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numXELE_A2+17+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile B #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4  ] } { incr k 1 } {

        remove element [expr int(($numXELE_A1+$i)*$numYele*$numZele+($numYele/2.0+$numXele11/2.0+4+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G #############
######################################################

set  numXELE_G    [expr $numXELE_d-1]

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4  ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G1 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+3+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G2 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+7+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G3 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+10+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G4 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+13+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}

######################################################
###### Removing Soil Elements for Pile G5 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4 ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numXELE_A2+17+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}


######################################################
###### Removing Soil Elements for Pile H #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_abut-$Df_abut/$zSize4  ] } { incr k 1 } {

        remove element [expr int(($numXELE_G+$i)*$numYele*$numZele+($numYele/2.0+$numXele11/2.0+4+$j)*$numZele+($HZ-$LP_abut)/$zSize1+$k)]
 
   }
  }
}


############################################################
########### Removing Soil Nodes for Pile A,B,G,H ###########
############################################################

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4-1] } { incr i 1 } {

     remove sp   [expr $p2_a+$numZnode+$i] 1
     remove sp   [expr $p2_a+$numZnode+$i] 2
     remove sp   [expr $p2_a+$numZnode+$i] 3
     remove node [expr $p2_a+$numZnode+$i]

	 remove sp   [expr $p2_a1+$numZnode+$i] 1
     remove sp   [expr $p2_a1+$numZnode+$i] 2
     remove sp   [expr $p2_a1+$numZnode+$i] 3
     remove node [expr $p2_a1+$numZnode+$i]
	 
	 remove sp   [expr $p2_a2+$numZnode+$i] 1
     remove sp   [expr $p2_a2+$numZnode+$i] 2
     remove sp   [expr $p2_a2+$numZnode+$i] 3
     remove node [expr $p2_a2+$numZnode+$i]
	 
	 remove sp   [expr $p2_a3+$numZnode+$i] 1
     remove sp   [expr $p2_a3+$numZnode+$i] 2
     remove sp   [expr $p2_a3+$numZnode+$i] 3
     remove node [expr $p2_a3+$numZnode+$i]
	 
	 remove sp   [expr $p2_a4+$numZnode+$i] 1
     remove sp   [expr $p2_a4+$numZnode+$i] 2
     remove sp   [expr $p2_a4+$numZnode+$i] 3
     remove node [expr $p2_a4+$numZnode+$i]
	 
	 remove sp   [expr $p2_a5+$numZnode+$i] 1
     remove sp   [expr $p2_a5+$numZnode+$i] 2
     remove sp   [expr $p2_a5+$numZnode+$i] 3
     remove node [expr $p2_a5+$numZnode+$i]
	 
     remove sp   [expr $p2_b+$numZnode+$i] 1
     remove sp   [expr $p2_b+$numZnode+$i] 2
     remove sp   [expr $p2_b+$numZnode+$i] 3
     remove node [expr $p2_b+$numZnode+$i]

     remove sp   [expr $p2_g+$numZnode+$i] 1
     remove sp   [expr $p2_g+$numZnode+$i] 2
     remove sp   [expr $p2_g+$numZnode+$i] 3
     remove node [expr $p2_g+$numZnode+$i]

	 remove sp   [expr $p2_g1+$numZnode+$i] 1
     remove sp   [expr $p2_g1+$numZnode+$i] 2
     remove sp   [expr $p2_g1+$numZnode+$i] 3
     remove node [expr $p2_g1+$numZnode+$i]
	 
	 remove sp   [expr $p2_g2+$numZnode+$i] 1
     remove sp   [expr $p2_g2+$numZnode+$i] 2
     remove sp   [expr $p2_g2+$numZnode+$i] 3
     remove node [expr $p2_g2+$numZnode+$i]
	 
	 remove sp   [expr $p2_g3+$numZnode+$i] 1
     remove sp   [expr $p2_g3+$numZnode+$i] 2
     remove sp   [expr $p2_g3+$numZnode+$i] 3
     remove node [expr $p2_g3+$numZnode+$i]
	 
	 remove sp   [expr $p2_g4+$numZnode+$i] 1
     remove sp   [expr $p2_g4+$numZnode+$i] 2
     remove sp   [expr $p2_g4+$numZnode+$i] 3
     remove node [expr $p2_g4+$numZnode+$i]
	 
	 remove sp   [expr $p2_g5+$numZnode+$i] 1
     remove sp   [expr $p2_g5+$numZnode+$i] 2
     remove sp   [expr $p2_g5+$numZnode+$i] 3
     remove node [expr $p2_g5+$numZnode+$i]
	 
     remove sp   [expr $p2_h+$numZnode+$i] 1
     remove sp   [expr $p2_h+$numZnode+$i] 2
     remove sp   [expr $p2_h+$numZnode+$i] 3
     remove node [expr $p2_h+$numZnode+$i]
    
    }

puts $meshFile "Elements"

########################################################
############### Pile Elements for PILE A ################
########################################################

set np_dum $np_a
set   rec_for_pile_A   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE A1 ################
########################################################

set np_dum $np_a1
set   rec_for_pile_A1   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE A2 ################
########################################################

set np_dum $np_a2
set   rec_for_pile_A2   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE A3 ################
########################################################

set np_dum $np_a3
set   rec_for_pile_A3   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE A4 ################
########################################################

set np_dum $np_a4
set   rec_for_pile_A4   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE A5 ################
########################################################

set np_dum $np_a5
set   rec_for_pile_A5   [expr $eleNum+1]
for { set k 1 } { $k <=  $pile_elenum_abut  } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
	
########################################################
############### Pile Elements for PILE B ################
########################################################

set np_dum $np_b
set   rec_for_pile_B   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }


########################################################
############### Pile Elements for PILE G ###############
########################################################

set np_dum $np_g
set   rec_for_pile_G   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE G1 ###############
########################################################

set np_dum $np_g1
set   rec_for_pile_G1   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################	
############### Pile Elements for PILE G2 ###############
########################################################

set np_dum $np_g2
set   rec_for_pile_G2   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE G3 ###############
########################################################

set np_dum $np_g3
set   rec_for_pile_G3   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE G4 ###############
########################################################

set np_dum $np_g4
set   rec_for_pile_G4   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE G5 ###############
########################################################

set np_dum $np_g5
set   rec_for_pile_G5   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts  1  101  -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE H ###############
########################################################

set np_dum $np_h
set   rec_for_pile_H   [expr $eleNum+1]

for { set k 1 } { $k <= $pile_elenum_abut } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101   -mass $massDens_pileabut
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }


# Define Horizontal Rigid Elements

###########################################################
########### Connecting Elements for PILE A ################
###########################################################

set ne6 $nodeNum149
set ne5 $np_a
set rec_connectionele_A  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

###########################################################
########### Connecting Elements for PILE A1 ################
###########################################################

set ne6 $nodeNum151
set ne5 $np_a1
set rec_connectionele_A1  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }
	
	
###########################################################
########### Connecting Elements for PILE A2 ################
###########################################################

set ne6 $nodeNum153
set ne5 $np_a2
set rec_connectionele_A2  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }
	
###########################################################
########### Connecting Elements for PILE A3 ################
###########################################################

set ne6 $nodeNum155
set ne5 $np_a3
set rec_connectionele_A3  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }	
	
	
###########################################################
########### Connecting Elements for PILE A4 ################
###########################################################

set ne6 $nodeNum157
set ne5 $np_a4
set rec_connectionele_A4  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }	
	
	
###########################################################
########### Connecting Elements for PILE A5 ################
###########################################################

set ne6 $nodeNum159
set ne5 $np_a5
set rec_connectionele_A5  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }	

	
###########################################################
########### Connecting Elements for PILE B ################
###########################################################

set ne6 $nodeNum161
set ne5 $np_b
set rec_connectionele_B  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }    
 
###########################################################
########### Connecting Elements for PILE G ################
###########################################################

set ne6      $nodeNum163
set ne5      $np_g
set rec_connectionele_G  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

###########################################################
########### Connecting Elements for PILE G1 ################
###########################################################

set ne6      $nodeNum165
set ne5      $np_g1
set rec_connectionele_G1  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE G2 ################
###########################################################

set ne6      $nodeNum167
set ne5      $np_g2
set rec_connectionele_G2  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE G3 ################
###########################################################

set ne6      $nodeNum169
set ne5      $np_g3
set rec_connectionele_G3  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE G4 ################
###########################################################

set ne6      $nodeNum171
set ne5      $np_g4
set rec_connectionele_G4  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE G5 ################
###########################################################

set ne6      $nodeNum173
set ne5      $np_g5
set rec_connectionele_G5  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

	
		
###########################################################
########### Connecting Elements for PILE H ################
###########################################################

set ne6      $nodeNum175
set ne5      $np_h
set rec_connectionele_H  $eleNum

for { set i 1 } { $i <=  [expr $pile_elenum_abut-$Df_abut/$zSize4 +1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*1)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*2)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_abut+1)*1+($pile_elenum_abut-$Df_abut/$zSize4+1)*3)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

puts $meshFile "end elements" 
  