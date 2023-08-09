
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates pile groups under bridge pier per site condition at Meloland Road Overpass, in El Centro, CA.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================


###############################################
###### Develop Pile Nodes at Bridge Pier ######
###############################################

# Creates input file for GiD software
puts $meshFile "MESH beamcolumn dimension 3 ElemType linear Nnode 2"
puts $meshFile "Coordinates"

set pile_elenum_pier  [expr int(($Czone3-$Df)/$zSize3+$Czone2/$zSize2+($LP_pier+$Df-$Czone3-$Czone2)/$zSize1)]

# See my Ph.D. dissertation for the pile labeling
########################################################
#################### PILE 1 ############################
########################################################
set   nodeNum7    $nodeNum6
set n $p1_1

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_1 }
  if { $i == 3 } { set n $p3_1 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum7 $xdim $ydim $zdim
      puts $meshFile "$nodeNum7 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_1 $nodeNum7       
      }
      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum7 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum7 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum7 1 2 3 
         }
      set nodeNum7 [expr $nodeNum7+1] 
    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum8    $nodeNum7
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum8 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
	  if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum8  [expr $np_1+$r-1]  1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum8  [expr $np_1+$r-1]  1 2 3 
	  }
      puts $meshFile "$nodeNum8 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum8  [expr $nodeNum8+1]
    } 

for { set i $np_1 } { $i <= [expr int($np_1+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE 2 ############################
########################################################
set   nodeNum101    $nodeNum8
set n $p1_2

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_2 }
  if { $i == 3 } { set n $p3_2 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum101 $xdim $ydim $zdim
      puts $meshFile "$nodeNum101 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_2 $nodeNum101       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum101 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum101 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum101 1 2 3 
         }

      set nodeNum101 [expr $nodeNum101+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum102    $nodeNum101

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum102 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum102  [expr $np_2+$r-1]  1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum102  [expr $np_2+$r-1]  1 2 3 
	  }
      puts $meshFile "$nodeNum102 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum102  [expr $nodeNum102+1]
    } 
	
for { set i $np_2 } { $i <= [expr int($np_2+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 3 ############################
########################################################
set   nodeNum103    $nodeNum102
set n $p1_3

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_3 }
  if { $i == 3 } { set n $p3_3 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum103 $xdim $ydim $zdim
      puts $meshFile "$nodeNum103 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_3 $nodeNum103       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum103 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum103 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum103 1 2 3 
         }

      set nodeNum103 [expr $nodeNum103+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum104   $nodeNum103
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum104 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum104  [expr $np_3+$r-1]   1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum104  [expr $np_3+$r-1]   1 2 3 
	  }
	  
 
      puts $meshFile "$nodeNum104 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum104  [expr $nodeNum104+1]
    } 

for { set i $np_3 } { $i <= [expr int($np_3+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 4 ############################
########################################################
set   nodeNum105    $nodeNum104
set n $p1_4

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_4 }
  if { $i == 3 } { set n $p3_4 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum105 $xdim $ydim $zdim
      puts $meshFile "$nodeNum105 $xdim $ydim $zdim"

      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_4 $nodeNum105       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum105 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum105 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum105 1 2 3 
         }

      set nodeNum105 [expr $nodeNum105+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum106   $nodeNum105

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum106 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ] $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum106  [expr $np_4+$r-1]  1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum106  [expr $np_4+$r-1]  1 2 3 
	  }

 
      puts $meshFile "$nodeNum106 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum106  [expr $nodeNum106+1]
    } 

for { set i $np_4 } { $i <= [expr int($np_4+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	

########################################################
#################### PILE 5 ############################
########################################################
set   nodeNum107    $nodeNum106
set n $p1_5

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_5 }
  if { $i == 3 } { set n $p3_5 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum107 $xdim $ydim $zdim
      puts $meshFile "$nodeNum107 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_5 $nodeNum107       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum107 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum107 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum107 1 2 3 
         }

      set nodeNum107 [expr $nodeNum107+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum108   $nodeNum107

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum108 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum108  [expr $np_5+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum108  [expr $np_5+$r-1]  1 2 3 
	  }
	   
 
      puts $meshFile "$nodeNum108 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0 ]  $zdim"
      set  nodeNum108  [expr $nodeNum108+1]
    } 

for { set i $np_5 } { $i <= [expr int($np_5+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	

########################################################
#################### PILE 6 ############################
########################################################
set   nodeNum109    $nodeNum108
set n $p1_6

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_6 }
  if { $i == 3 } { set n $p3_6 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum109 $xdim $ydim $zdim
      puts $meshFile "$nodeNum109 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_6 $nodeNum109       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum109  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum109  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum109  1 2 3 
         }

      set nodeNum109 [expr $nodeNum109+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum110    $nodeNum109

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum110 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum110  [expr $np_6+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum110  [expr $np_6+$r-1]  1 2 3 
	  }

      puts $meshFile "$nodeNum110  [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum110  [expr $nodeNum110+1]
    } 

for { set i $np_6 } { $i <= [expr int($np_6+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE 7 ############################
########################################################
set   nodeNum111    $nodeNum110
set n $p1_7

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_7 }
  if { $i == 3 } { set n $p3_7 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum111 $xdim $ydim $zdim
      puts $meshFile "$nodeNum111 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_7 $nodeNum111       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum111 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum111 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum111 1 2 3 
         }

      set nodeNum111 [expr $nodeNum111+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum112    $nodeNum111

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum112 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim

      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum112  [expr $np_7+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum112  [expr $np_7+$r-1]  1 2 3 
	  }
	   
      puts $meshFile "$nodeNum112 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum112  [expr $nodeNum112+1]
    } 

for { set i $np_7 } { $i <= [expr int($np_7+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 8 ############################
########################################################
set   nodeNum113    $nodeNum112
set n $p1_8

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_8 }
  if { $i == 3 } { set n $p3_8 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum113 $xdim $ydim $zdim
      puts $meshFile "$nodeNum113 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_8 $nodeNum113       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum113 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum113 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum113 1 2 3 
         }

      set nodeNum113 [expr $nodeNum113+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum114   $nodeNum113

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum114 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum114  [expr $np_8+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum114  [expr $np_8+$r-1] 1 2 3 
	  }
	
      puts $meshFile "$nodeNum114 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum114  [expr $nodeNum114+1]
    } 

for { set i $np_8 } { $i <= [expr int($np_8+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 9 ############################
########################################################
set   nodeNum115    $nodeNum114
set n $p1_9

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_9 }
  if { $i == 3 } { set n $p3_9 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum115 $xdim $ydim $zdim
      puts $meshFile "$nodeNum115 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_9 $nodeNum115       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum115 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum115 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum115 1 2 3 
         }

      set nodeNum115 [expr $nodeNum115+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum116   $nodeNum115

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum116 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ] $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum116  [expr $np_9+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum116  [expr $np_9+$r-1] 1 2 3 
	  }
	
      puts $meshFile "$nodeNum116 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum116  [expr $nodeNum116+1]
    } 

for { set i $np_9 } { $i <= [expr int($np_9+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	
	
########################################################
#################### PILE 10 ############################
########################################################
set   nodeNum117    $nodeNum116
set n $p1_10

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_10 }
  if { $i == 3 } { set n $p3_10 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
     set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$WInfemb+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum117 $xdim $ydim $zdim
      puts $meshFile "$nodeNum117 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_10 $nodeNum117       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum117 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum117 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum117 1 2 3 
         }

      set nodeNum117 [expr $nodeNum117+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum118   $nodeNum117

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum118 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum118  [expr $np_10+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum118  [expr $np_10+$r-1] 1 2 3 
	  }
	 
      puts $meshFile "$nodeNum118 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0 ]  $zdim"
      set  nodeNum118  [expr $nodeNum118+1]
    } 

for { set i $np_10 } { $i <= [expr int($np_10+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	
		
########################################################
#################### PILE 11 ############################
########################################################
set   nodeNum119    $nodeNum118
set n $p1_11

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_11 }
  if { $i == 3 } { set n $p3_11 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum119 $xdim $ydim $zdim
      puts $meshFile "$nodeNum119 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_11 $nodeNum119       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum119  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum119  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum119  1 2 3 
         }

      set nodeNum119 [expr $nodeNum119+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum120    $nodeNum119

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
	  
      node  $nodeNum120 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum120  [expr $np_11+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum120  [expr $np_11+$r-1] 1 2 3 
	  }
	   
      puts $meshFile "$nodeNum120  [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum120  [expr $nodeNum120+1]
    } 

for { set i $np_11 } { $i <= [expr int($np_11+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE 12 ############################
########################################################
set   nodeNum121    $nodeNum120
set n $p1_12

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_12 }
  if { $i == 3 } { set n $p3_12 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
  
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }     
      node $nodeNum121 $xdim $ydim $zdim
      puts $meshFile "$nodeNum121 $xdim $ydim $zdim"
      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_12 $nodeNum121       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum121 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum121 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum121 1 2 3 
         }

      set nodeNum121 [expr $nodeNum121+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum122    $nodeNum121

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
                                                                                                              	  
      node  $nodeNum122 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum122  [expr $np_12+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum122  [expr $np_12+$r-1] 1 2 3 
	  }
	  
 
      puts $meshFile "$nodeNum122 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum122  [expr $nodeNum122+1]
    } 

for { set i $np_12 } { $i <= [expr int($np_12+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }
	
########################################################
#################### PILE 13 ############################
########################################################
set   nodeNum123    $nodeNum122
set n $p1_13

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_13 }
  if { $i == 3 } { set n $p3_13 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum123 $xdim $ydim $zdim
      puts $meshFile "$nodeNum123 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_13 $nodeNum123       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum123 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum123 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum123 1 2 3 
         }

      set nodeNum123 [expr $nodeNum123+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum124   $nodeNum123
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum124 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
     equalDOF  $nodeNum124  [expr $np_13+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum124  [expr $np_13+$r-1] 1 2 3 
	  }
	 
      puts $meshFile "$nodeNum124 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum124  [expr $nodeNum124+1]
    } 

for { set i $np_13 } { $i <= [expr int($np_13+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 14 ############################
########################################################
set   nodeNum125    $nodeNum124
set n $p1_14

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_14 }
  if { $i == 3 } { set n $p3_14 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum125 $xdim $ydim $zdim
      puts $meshFile "$nodeNum125 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_14 $nodeNum125       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum125 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum125 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum125 1 2 3 
         }

      set nodeNum125 [expr $nodeNum125+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum126   $nodeNum125

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum126 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ] $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
     equalDOF  $nodeNum126  [expr $np_14+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum126  [expr $np_14+$r-1] 1 2 3 
	  }
	  
 
      puts $meshFile "$nodeNum126 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum126  [expr $nodeNum126+1]
    } 

for { set i $np_14 } { $i <= [expr int($np_14+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	
	
########################################################
#################### PILE 15 ############################
########################################################
set   nodeNum127    $nodeNum126
set n $p1_15

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_15 }
  if { $i == 3 } { set n $p3_15 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
     set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum127 $xdim $ydim $zdim
      puts $meshFile "$nodeNum127 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_15 $nodeNum127       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum127 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum127 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum127 1 2 3 
         }

      set nodeNum127 [expr $nodeNum127+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum128   $nodeNum127
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum128 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum128  [expr $np_15+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum128  [expr $np_15+$r-1] 1 2 3 
	  }
	 
      puts $meshFile "$nodeNum128 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0 ]  $zdim"
      set  nodeNum128  [expr $nodeNum128+1]
    } 

for { set i $np_15 } { $i <= [expr int($np_15+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	

########################################################
#################### PILE 16 ############################
########################################################
set   nodeNum129    $nodeNum128
set n $p1_16

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_16 }
  if { $i == 3 } { set n $p3_16 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum129 $xdim $ydim $zdim
      puts $meshFile "$nodeNum129 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_16 $nodeNum129       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum129  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum129  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum129  1 2 3 
         }

      set nodeNum129 [expr $nodeNum129+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum130    $nodeNum129

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum130 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
     if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum130  [expr $np_16+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum130  [expr $np_16+$r-1] 1 2 3 
	  }
      puts $meshFile "$nodeNum130  [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum130  [expr $nodeNum130+1]
    } 

for { set i $np_16 } { $i <= [expr int($np_16+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 


########################################################
#################### PILE 17 ############################
########################################################
set   nodeNum131    $nodeNum130
set n $p1_17

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_17 }
  if { $i == 3 } { set n $p3_17 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum131 $xdim $ydim $zdim
      puts $meshFile "$nodeNum131 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_17 $nodeNum131      
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum131 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum131 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum131 1 2 3 
         }

      set nodeNum131 [expr $nodeNum131+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum132    $nodeNum131

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum132 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      
	  if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum132  [expr $np_17+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum132  [expr $np_17+$r-1] 1 2 3 
	  }
 
      puts $meshFile "$nodeNum132 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum132  [expr $nodeNum132+1]
    } 

for { set i $np_17 } { $i <= [expr int($np_17+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }
	
########################################################
#################### PILE 18 ############################
########################################################
set   nodeNum133    $nodeNum132
set n $p1_18

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_18 }
  if { $i == 3 } { set n $p3_18 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum133 $xdim $ydim $zdim
      puts $meshFile "$nodeNum133 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_18 $nodeNum133       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum133 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum133 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum133 1 2 3 
         }

      set nodeNum133 [expr $nodeNum133+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum134   $nodeNum133

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum134 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum134  [expr $np_18+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum134  [expr $np_18+$r-1] 1 2 3 
	  }
 
      puts $meshFile "$nodeNum134 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum134  [expr $nodeNum134+1]
    } 

for { set i $np_18 } { $i <= [expr int($np_18+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 19 ############################
########################################################
set   nodeNum135    $nodeNum134
set n $p1_19

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_19 }
  if { $i == 3 } { set n $p3_19 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum135 $xdim $ydim $zdim
      puts $meshFile "$nodeNum135 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_19 $nodeNum135       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum135 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum135 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum135 1 2 3 
         }

      set nodeNum135 [expr $nodeNum135+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum136   $nodeNum135

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum136 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ] $zdim

      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum136  [expr $np_19+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum136  [expr $np_19+$r-1] 1 2 3 
	  } 
 
      puts $meshFile "$nodeNum136 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum136  [expr $nodeNum136+1]
    } 

for { set i $np_19 } { $i <= [expr int($np_19+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	
	
########################################################
#################### PILE 20 ############################
########################################################
set   nodeNum137    $nodeNum136
set n $p1_20

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_20 }
  if { $i == 3 } { set n $p3_20 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum137 $xdim $ydim $zdim
      puts $meshFile "$nodeNum137 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_20 $nodeNum137       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum137 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum137 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum137 1 2 3 
         }

      set nodeNum137 [expr $nodeNum137+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum138   $nodeNum137
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum138 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      
	  if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum138  [expr $np_20+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum138  [expr $np_20+$r-1] 1 2 3 
	  } 
	  
 
      puts $meshFile "$nodeNum138 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0 ]  $zdim"
      set  nodeNum138  [expr $nodeNum138+1]
    } 

for { set i $np_20 } { $i <= [expr int($np_20+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	

########################################################
#################### PILE 21 ############################
########################################################
set   nodeNum139    $nodeNum138
set n $p1_21

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_21 }
  if { $i == 3 } { set n $p3_21 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum139 $xdim $ydim $zdim
      puts $meshFile "$nodeNum139 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_21 $nodeNum139       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum139  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum139  1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum139  1 2 3 
         }

      set nodeNum139 [expr $nodeNum139+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum140    $nodeNum139

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum140 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum140  [expr $np_21+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum140  [expr $np_21+$r-1] 1 2 3 
	  } 	 
 
      puts $meshFile "$nodeNum140  [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile/2.0 ]  $zdim"
      set  nodeNum140  [expr $nodeNum140+1]
    } 

for { set i $np_21 } { $i <= [expr int($np_21+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    } 

########################################################
#################### PILE 22 ############################
########################################################
set   nodeNum141    $nodeNum140

set n $p1_22

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_22 }
  if { $i == 3 } { set n $p3_22 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum141 $xdim $ydim $zdim
      puts $meshFile "$nodeNum141 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_22 $nodeNum141      
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum141 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum141 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum141 1 2 3 
         }

      set nodeNum141 [expr $nodeNum141+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum142    $nodeNum141
for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum142 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim
      
	  if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum142  [expr $np_22+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum142  [expr $np_22+$r-1] 1 2 3 
	  } 
 
      puts $meshFile "$nodeNum142 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum142  [expr $nodeNum142+1]
    } 

for { set i $np_22 } { $i <= [expr int($np_22+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }
	
########################################################
#################### PILE 23 ############################
########################################################
set   nodeNum143    $nodeNum142
set n $p1_23

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_23 }
  if { $i == 3 } { set n $p3_23 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
     set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($i-1)*$xSize10 ]
     set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum143 $xdim $ydim $zdim
      puts $meshFile "$nodeNum143 $xdim $ydim $zdim"
    
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_23 $nodeNum143       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum143 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum143 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum143 1 2 3 
         }

      set nodeNum143 [expr $nodeNum143+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum144   $nodeNum143

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum144 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum144  [expr $np_23+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum144  [expr $np_23+$r-1] 1 2 3 
	  } 
	   
 
      puts $meshFile "$nodeNum144 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum144  [expr $nodeNum144+1]
    } 

for { set i $np_23 } { $i <= [expr int($np_23+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }

########################################################
#################### PILE 24 ############################
########################################################
set   nodeNum145    $nodeNum144

set n $p1_24

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_24 }
  if { $i == 3 } { set n $p3_24 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum145 $xdim $ydim $zdim
      puts $meshFile "$nodeNum145 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_24 $nodeNum145       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum145 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum145 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum145 1 2 3 
         }

      set nodeNum145 [expr $nodeNum145+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum146   $nodeNum145

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }

      node  $nodeNum146 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ] $zdim
      
	  if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum146  [expr $np_24+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum146  [expr $np_24+$r-1] 1 2 3 
	  } 
 
      puts $meshFile "$nodeNum146 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0 ]  $zdim"
      set  nodeNum146  [expr $nodeNum146+1]
    } 

for { set i $np_24 } { $i <= [expr int($np_24+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }	
	
########################################################
#################### PILE 25 ############################
########################################################
set   nodeNum147    $nodeNum146
set n $p1_25

for {set i 1} {$i <= 3} {incr i 1} {
  if { $i == 2 } { set n $p2_25 }
  if { $i == 3 } { set n $p3_25 }
  for {set j 1} {$j <= 3} {incr j 1} {
    for {set k 1} {$k <= [expr $pile_elenum_pier+1]} {incr k 1} {
 
      set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($i-1)*$xSize10 ]
      set ydim [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$WInfemb+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+($j-1)*$xSize10 ]

      if { $k <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($k-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $k > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($k-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      
      node $nodeNum147 $xdim $ydim $zdim
      puts $meshFile "$nodeNum147 $xdim $ydim $zdim"

      
      if { $k == 1 && $i == 2 && $j == 2 } {
      set np_25 $nodeNum147       
      }

      #######  Slaving Pile Nodes to Soil Nodes
      if { $i != 2 || $j != 2 } {
          equalDOF  [expr $n+$k-1] $nodeNum147 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == 1 } {
          equalDOF  [expr $n+$k-1] $nodeNum147 1 2 3 
         } 
      if { $i == 2 && $j == 2 && $k == [expr $pile_elenum_pier+1] } {
          equalDOF  [expr $n+$k-1] $nodeNum147 1 2 3 
         }

      set nodeNum147 [expr $nodeNum147+1] 

    }     
    set n [expr ($n+$numZnode)]
  }
}

set nodeNum148   $nodeNum147

for { set r 1 } { $r <= [expr ($pile_elenum_pier+1)] } { incr r 1 } {

      if { $r <= [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1] } {
          set   zdim [expr ($r-1)*$zSize1 + $HZ   -   ($LP_pier+$Df)]
      } else {
        if { $r > [expr ($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1] } {
          set   zdim [expr $Czone1+$Czone2+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+$numZele2+1))*$zSize3]
      } else {
          set   zdim [expr $Czone1+($r-(($LP_pier+$Df-$Czone3-$Czone2)/$zSize1+1))*$zSize2]
      }
      }
      node  $nodeNum148 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0  ]  $zdim
      if { $r == [expr ($pile_elenum_pier+1)] } {
      equalDOF  $nodeNum148  [expr $np_25+$r-1] 1 2 3 4 5 6
      } else {
	  equalDOF  $nodeNum148  [expr $np_25+$r-1] 1 2 3 
	  } 
 
      puts $meshFile "$nodeNum148 [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0] [expr $Bzone1+$Bzone2+$Bzone3+2.0*$WInfemb+$BP_emb+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$WInfpile+$BP_pile/2.0 ]  $zdim"
      set  nodeNum148  [expr $nodeNum148+1]
    } 

for { set i $np_25 } { $i <= [expr int($np_25+$pile_elenum_pier)]} { incr i 1 } {
    fix $i 0 0 0 0 0 1
    }		
	
puts $meshFile "end coordinates"

######################################################
###### Removing Soil Elements for Pile 1 #############
######################################################

set  numXELE_C    [expr $numELE_b+$numXele8+$numXele9+1]
set  numXELE_D    [expr $numYele1 +$numYele2 +$numYele3+3]

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
        ####puts "remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]"		
 }
  }
}

######################################################
###### Removing Soil Elements for Pile 2 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+5+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 3 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+8+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 4 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+11+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 5 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+$i)*$numYele*$numZele+($numXELE_D+15+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 6 #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+4+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 7 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+4+$i)*$numYele*$numZele+($numXELE_D+5+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 8 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+4+$i)*$numYele*$numZele+($numXELE_D+8+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 9 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+4+$i)*$numYele*$numZele+($numXELE_D+11+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 10 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+4+$i)*$numYele*$numZele+($numXELE_D+15+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 11 #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+7+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 12 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+7+$i)*$numYele*$numZele+($numXELE_D+5+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 13 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+7+$i)*$numYele*$numZele+($numXELE_D+8+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 14 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+7+$i)*$numYele*$numZele+($numXELE_D+11+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 15 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+7+$i)*$numYele*$numZele+($numXELE_D+15+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 16 #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+10+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 17 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+10+$i)*$numYele*$numZele+($numXELE_D+5+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 18 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+10+$i)*$numYele*$numZele+($numXELE_D+8+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 19 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+10+$i)*$numYele*$numZele+($numXELE_D+11+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 20 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+10+$i)*$numYele*$numZele+($numXELE_D+15+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}


######################################################
###### Removing Soil Elements for Pile 21 #############
######################################################


for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+14+$i)*$numYele*$numZele+($numXELE_D+1+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 22 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+14+$i)*$numYele*$numZele+($numXELE_D+5+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 23 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+14+$i)*$numYele*$numZele+($numXELE_D+8+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 24 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+14+$i)*$numYele*$numZele+($numXELE_D+11+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

######################################################
###### Removing Soil Elements for Pile 25 #############
######################################################

for { set i 0 } { $i <= 1 } { incr i 1 } {
  for { set j 0 } { $j <= 1 } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $pile_elenum_pier ] } { incr k 1 } {

        remove element [expr int(($numXELE_C+14+$i)*$numYele*$numZele+($numXELE_D+15+$j)*$numZele+($HZ-($LP_pier+$Df))/$zSize1+$k)]
   }
  }
}

###########################################################
########## Removing Soil Nodes for Piles 1 to 25 ###########
###########################################################

for { set i 1 } { $i <=  [expr $pile_elenum_pier-1] } { incr i 1 } {

     remove sp   [expr $p2_1+$numZnode+$i] 1
     remove sp   [expr $p2_1+$numZnode+$i] 2
     remove sp   [expr $p2_1+$numZnode+$i] 3
     remove node [expr $p2_1+$numZnode+$i]
     
     remove sp   [expr $p2_2+$numZnode+$i] 1
     remove sp   [expr $p2_2+$numZnode+$i] 2
     remove sp   [expr $p2_2+$numZnode+$i] 3
     remove node [expr $p2_2+$numZnode+$i]

     remove sp   [expr $p2_3+$numZnode+$i] 1
     remove sp   [expr $p2_3+$numZnode+$i] 2
     remove sp   [expr $p2_3+$numZnode+$i] 3
     remove node [expr $p2_3+$numZnode+$i]

     remove sp   [expr $p2_4+$numZnode+$i] 1
     remove sp   [expr $p2_4+$numZnode+$i] 2
     remove sp   [expr $p2_4+$numZnode+$i] 3
     remove node [expr $p2_4+$numZnode+$i]
	 
	 remove sp   [expr $p2_5+$numZnode+$i] 1
     remove sp   [expr $p2_5+$numZnode+$i] 2
     remove sp   [expr $p2_5+$numZnode+$i] 3
     remove node [expr $p2_5+$numZnode+$i]
	 
 
     remove sp   [expr $p2_6+$numZnode+$i] 1
     remove sp   [expr $p2_6+$numZnode+$i] 2
     remove sp   [expr $p2_6+$numZnode+$i] 3
     remove node [expr $p2_6+$numZnode+$i]

     remove sp   [expr $p2_7+$numZnode+$i] 1
     remove sp   [expr $p2_7+$numZnode+$i] 2
     remove sp   [expr $p2_7+$numZnode+$i] 3
     remove node [expr $p2_7+$numZnode+$i]

     remove sp   [expr $p2_8+$numZnode+$i] 1
     remove sp   [expr $p2_8+$numZnode+$i] 2
     remove sp   [expr $p2_8+$numZnode+$i] 3
     remove node [expr $p2_8+$numZnode+$i]

     remove sp   [expr $p2_9+$numZnode+$i] 1
     remove sp   [expr $p2_9+$numZnode+$i] 2
     remove sp   [expr $p2_9+$numZnode+$i] 3
     remove node [expr $p2_9+$numZnode+$i]
	 
	 remove sp   [expr $p2_10+$numZnode+$i] 1
     remove sp   [expr $p2_10+$numZnode+$i] 2
     remove sp   [expr $p2_10+$numZnode+$i] 3
     remove node [expr $p2_10+$numZnode+$i]

     remove sp   [expr $p2_11+$numZnode+$i] 1
     remove sp   [expr $p2_11+$numZnode+$i] 2
     remove sp   [expr $p2_11+$numZnode+$i] 3
     remove node [expr $p2_11+$numZnode+$i]

     remove sp   [expr $p2_12+$numZnode+$i] 1
     remove sp   [expr $p2_12+$numZnode+$i] 2
     remove sp   [expr $p2_12+$numZnode+$i] 3
     remove node [expr $p2_12+$numZnode+$i]

     remove sp   [expr $p2_13+$numZnode+$i] 1
     remove sp   [expr $p2_13+$numZnode+$i] 2
     remove sp   [expr $p2_13+$numZnode+$i] 3
     remove node [expr $p2_13+$numZnode+$i]

     remove sp   [expr $p2_14+$numZnode+$i] 1
     remove sp   [expr $p2_14+$numZnode+$i] 2
     remove sp   [expr $p2_14+$numZnode+$i] 3
     remove node [expr $p2_14+$numZnode+$i]
	 
	 remove sp   [expr $p2_15+$numZnode+$i] 1
     remove sp   [expr $p2_15+$numZnode+$i] 2
     remove sp   [expr $p2_15+$numZnode+$i] 3
     remove node [expr $p2_15+$numZnode+$i]

     remove sp   [expr $p2_16+$numZnode+$i] 1
     remove sp   [expr $p2_16+$numZnode+$i] 2
     remove sp   [expr $p2_16+$numZnode+$i] 3
     remove node [expr $p2_16+$numZnode+$i]
	 
     remove sp   [expr $p2_17+$numZnode+$i] 1
     remove sp   [expr $p2_17+$numZnode+$i] 2
     remove sp   [expr $p2_17+$numZnode+$i] 3
     remove node [expr $p2_17+$numZnode+$i]

     remove sp   [expr $p2_18+$numZnode+$i] 1
     remove sp   [expr $p2_18+$numZnode+$i] 2
     remove sp   [expr $p2_18+$numZnode+$i] 3
     remove node [expr $p2_18+$numZnode+$i]

     remove sp   [expr $p2_19+$numZnode+$i] 1
     remove sp   [expr $p2_19+$numZnode+$i] 2
     remove sp   [expr $p2_19+$numZnode+$i] 3
     remove node [expr $p2_19+$numZnode+$i]
	 
	 remove sp   [expr $p2_20+$numZnode+$i] 1
     remove sp   [expr $p2_20+$numZnode+$i] 2
     remove sp   [expr $p2_20+$numZnode+$i] 3
     remove node [expr $p2_20+$numZnode+$i]

	 remove sp   [expr $p2_21+$numZnode+$i] 1
     remove sp   [expr $p2_21+$numZnode+$i] 2
     remove sp   [expr $p2_21+$numZnode+$i] 3
     remove node [expr $p2_21+$numZnode+$i]
	 
     remove sp   [expr $p2_22+$numZnode+$i] 1
     remove sp   [expr $p2_22+$numZnode+$i] 2
     remove sp   [expr $p2_22+$numZnode+$i] 3
     remove node [expr $p2_22+$numZnode+$i]

     remove sp   [expr $p2_23+$numZnode+$i] 1
     remove sp   [expr $p2_23+$numZnode+$i] 2
     remove sp   [expr $p2_23+$numZnode+$i] 3
     remove node [expr $p2_23+$numZnode+$i]

     remove sp   [expr $p2_24+$numZnode+$i] 1
     remove sp   [expr $p2_24+$numZnode+$i] 2
     remove sp   [expr $p2_24+$numZnode+$i] 3
     remove node [expr $p2_24+$numZnode+$i]
	 
	 remove sp   [expr $p2_25+$numZnode+$i] 1
     remove sp   [expr $p2_25+$numZnode+$i] 2
     remove sp   [expr $p2_25+$numZnode+$i] 3
     remove node [expr $p2_25+$numZnode+$i]
	 
	 
   }

puts $meshFile "Elements"

geomTransf PDelta 101 1 0 0
geomTransf Linear 102 0 0 -1

set numIntgrPts 5

########################################################
############### Pile Elements for PILE 1 ###############
########################################################

set np_dum $np_1
set rec_for_pile_1 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }


########################################################
############### Pile Elements for PILE 2 ###############
########################################################

set np_dum $np_2
set rec_for_pile_2 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 3 ###############
########################################################

set np_dum $np_3
set rec_for_pile_3 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE 4 ###############
########################################################

set np_dum $np_4
set rec_for_pile_4 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	
	
########################################################
############### Pile Elements for PILE 5 ###############
########################################################

set np_dum $np_5
set rec_for_pile_5 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }		
	
########################################################
############### Pile Elements for PILE 6 ###############
########################################################

set np_dum $np_6
set rec_for_pile_6 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	
	
########################################################
############### Pile Elements for PILE 7 ###############
########################################################

set np_dum $np_7
set rec_for_pile_7 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }		
	
########################################################
############### Pile Elements for PILE 8 ###############
########################################################

set np_dum $np_8
set rec_for_pile_8 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE 9 ###############
########################################################

set np_dum $np_9
set rec_for_pile_9 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	

########################################################
############### Pile Elements for PILE 10 ################
########################################################

set np_dum $np_10
set rec_for_pile_10 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 11 ################
########################################################

set np_dum $np_11
set rec_for_pile_11 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 12 ################
########################################################

set np_dum $np_12
set rec_for_pile_12 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 13 ################
########################################################

set np_dum $np_13
set rec_for_pile_13 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 14 ################
########################################################

set np_dum $np_14
set rec_for_pile_14 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 15 ################
########################################################

set np_dum $np_15
set rec_for_pile_15 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 16 ################
########################################################

set np_dum $np_16
set rec_for_pile_16 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 17 ################
########################################################

set np_dum $np_17
set rec_for_pile_17 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 18 ################
########################################################

set np_dum $np_18
set rec_for_pile_18 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 19 ################
########################################################

set np_dum $np_19
set rec_for_pile_19 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 20 ################
########################################################

set np_dum $np_20
set rec_for_pile_20 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }

########################################################
############### Pile Elements for PILE 21 ################
########################################################

set np_dum $np_21
set rec_for_pile_21 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }
	
########################################################
############### Pile Elements for PILE 22 ################
########################################################

set np_dum $np_22
set rec_for_pile_22 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	
	
########################################################
############### Pile Elements for PILE 23 ################
########################################################

set np_dum $np_23
set rec_for_pile_23 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	

########################################################
############### Pile Elements for PILE 24 ################
########################################################

set np_dum $np_24
set rec_for_pile_24 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	

########################################################
############### Pile Elements for PILE 25 ################
########################################################

set np_dum $np_25
set rec_for_pile_25 [expr $eleNum+1]
for { set k 1 } { $k <= $pile_elenum_pier } { incr k 1 } { 

  set eleNum [expr $eleNum+1]
  element dispBeamColumn  $eleNum $np_dum [expr $np_dum+1]   $numIntgrPts 1  101 -mass $massDens_pilepier 
  puts $meshFile "$eleNum $np_dum [expr $np_dum+1] 10"
  set  np_dum   [expr $np_dum+1]

    }	


# Define Horizontal Rigid Connections
 
###########################################################
########### Connecting Elements for PILE 1 ################
###########################################################

set ne6 $nodeNum7
set ne5 $np_1
set rec_connectionele_1 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

###########################################################
########### Connecting Elements for PILE 2 ################
###########################################################

set ne6 $nodeNum101
set ne5 $np_2
set rec_connectionele_2 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  	

###########################################################
########### Connecting Elements for PILE 3 ################
###########################################################

set ne6 $nodeNum103
set ne5 $np_3
set rec_connectionele_3 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  	

###########################################################
########### Connecting Elements for PILE 4 ################
###########################################################

set ne6 $nodeNum105
set ne5 $np_4
set rec_connectionele_4 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  	

###########################################################
########### Connecting Elements for PILE 5 ################
###########################################################

set ne6 $nodeNum107
set ne5 $np_5
set rec_connectionele_5 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  	


###########################################################
########### Connecting Elements for PILE 6 ################
###########################################################

set ne6 $nodeNum109
set ne5 $np_6
set rec_connectionele_6 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

	
###########################################################
########### Connecting Elements for PILE 7 ################
###########################################################

set ne6 $nodeNum111
set ne5 $np_7
set rec_connectionele_7 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }  

###########################################################
########### Connecting Elements for PILE 8 ################
###########################################################

set ne6 $nodeNum113
set ne5 $np_8
set rec_connectionele_8 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE 9 ################
###########################################################

set ne6 $nodeNum115
set ne5 $np_9
set rec_connectionele_9 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 


###########################################################
########### Connecting Elements for PILE 10 ################
###########################################################

set ne6 $nodeNum117
set ne5 $np_10
set rec_connectionele_10 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 
	
###########################################################
########### Connecting Elements for PILE 11 ################
###########################################################

set ne6 $nodeNum119
set ne5 $np_11
set rec_connectionele_11 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 
	
###########################################################
########### Connecting Elements for PILE 12 ################
###########################################################

set ne6 $nodeNum121
set ne5 $np_12
set rec_connectionele_12 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 
	
###########################################################
########### Connecting Elements for PILE 13 ################
###########################################################

set ne6 $nodeNum123
set ne5 $np_13
set rec_connectionele_13 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE 14 ################
###########################################################

set ne6 $nodeNum125
set ne5 $np_14
set rec_connectionele_14 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    } 

###########################################################
########### Connecting Elements for PILE 15 ################
###########################################################

set ne6 $nodeNum127
set ne5 $np_15
set rec_connectionele_15 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 16 ################
###########################################################

set ne6 $nodeNum129
set ne5 $np_16
set rec_connectionele_16 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 17 ################
###########################################################

set ne6 $nodeNum131
set ne5 $np_17
set rec_connectionele_17 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 18 ################
###########################################################

set ne6 $nodeNum133
set ne5 $np_18
set rec_connectionele_18 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 19 ################
###########################################################

set ne6 $nodeNum135
set ne5 $np_19
set rec_connectionele_19 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 20 ################
###########################################################

set ne6 $nodeNum137
set ne5 $np_20
set rec_connectionele_20 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 21 ################
###########################################################

set ne6 $nodeNum139
set ne5 $np_21
set rec_connectionele_21 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }


###########################################################
########### Connecting Elements for PILE 22 ################
###########################################################

set ne6 $nodeNum141
set ne5 $np_22
set rec_connectionele_22 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 23 ################
###########################################################

set ne6 $nodeNum143
set ne5 $np_23
set rec_connectionele_23 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 24 ################
###########################################################

set ne6 $nodeNum145
set ne5 $np_24
set rec_connectionele_24 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }

###########################################################
########### Connecting Elements for PILE 25 ################
###########################################################

set ne6 $nodeNum147
set ne5 $np_25
set rec_connectionele_25 $eleNum
for { set i 1 } { $i <=  [expr $pile_elenum_pier+1]  } { incr i 1 } {

    element elasticBeamColumn    [expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102              
    element elasticBeamColumn    [expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102
    element elasticBeamColumn    [expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)]   $Acon $Econ $Gcon $Jcon $Iycon $Izcon 102

    puts $meshFile "[expr $eleNum+1]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*4)] 10"
    puts $meshFile "[expr $eleNum+2]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+3]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+4]   $ne6  [expr int($ne5 - ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+5]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*1)] 10"
    puts $meshFile "[expr $eleNum+6]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*2)] 10"
    puts $meshFile "[expr $eleNum+7]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*3)] 10"
    puts $meshFile "[expr $eleNum+8]   $ne6  [expr int($ne5 + ($pile_elenum_pier+1)*4)] 10"

    set ne6      [expr $ne6+1]
    set ne5      [expr $ne5+1]
    set eleNum   [expr $eleNum+8]

    }
	
puts $meshFile "end elements"
  
