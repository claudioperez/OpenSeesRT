
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates the bridge abutment backwall and wing walls per site condition at Meloland Road Overpass, in El Centro, CA.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================
 
##################################################
###### Develop Left Abutment Backwall Nodes ######
##################################################

# Creates input file for GiD software
puts $meshFile "MESH shell dimension 3 ElemType Quadrilateral Nnode 4"
puts $meshFile "Coordinates"

set  nodeNum26  $nodeNum25

for { set j 1 } { $j <= [ expr $numYeleabut+1] } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
    
         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0]

          if { $j <= [expr $numYele3 +1] } {
            set ydim [expr $Bzone1+$Bzone2+($j-1)*$ySize3]
            } else {
         if { $j == [expr $numYele3+1 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
            } else {
         if { $j == [expr $numYele3+2 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb/2.0] 
            } else {
         if { $j == [expr $numYele3+3 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb]                 
            } else {
         if { $j == [expr $numYele3+4 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
         if { $j == [expr $numYele3+5 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
             } else {
         if { $j == [expr $numYele3+6 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]                       
             } else {
         if { $j == [expr $numYele3+7 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]
             } else {
			 
         if { $j == [expr $numYele3+8+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+9+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+10+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]            
             } else {
         if { $j == [expr $numYele3+11+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+12+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+13+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]            
             } else {
         if { $j == [expr $numYele3+14+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+15+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+16+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]            
             } else {			 
         if { $j == [expr $numYele3+17+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]            
             } else {				 
         if { $j == [expr $numYele3+7+$numXele11+1+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]                                           
             } else {
         if { $j == [expr $numYele3+7+$numXele11+2+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]                                                                          
              } else {
         if { $j == [expr $numYele3+7+$numXele11+3+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]                                                                                                      
              } else {
         if { $j == [expr $numYele3+7+$numXele11+4+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]                                                                                                                                    
              } else {
         if { $j == [expr $numYele3+7+$numXele11+5+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]                                                                                                                                                                  
               } else {
         if { $j == [expr $numYele3+7+$numXele11+6+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
               } else {
         if { $j == [expr $numYele3+7+$numXele11+7+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]            
               } else {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb+($j-($numYele3+7+$numXele11+7+ 1))*$ySize5] 
               }}}}}}}}}}}}}}}}}}}}}}}}}

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum26 $xdim $ydim $zdim
         puts $meshFile "$nodeNum26 $xdim $ydim $zdim"
 
         fix $nodeNum26 0 0 0 0 0 1
 
         if { $k == 1 && $j == [expr int($numYele3+1+1)] } {
            equalDOF [expr int($np_a+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+5+1)] } {
            equalDOF [expr int($np_a1+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+9+1)] } {
            equalDOF [expr int($np_a2+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+12+1)] } {
            equalDOF [expr int($np_a3+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+15+1)] } {
            equalDOF [expr int($np_a4+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+19+1)] } {
            equalDOF [expr int($np_a5+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5 
            } else {			
         if { $k == 1 && $j == [expr int($numYele3+7+$numXele11+6+1)] } {
            equalDOF [expr int($np_b+$pile_elenum_abut)] $nodeNum26 1 2 3 4 5
            } else {
           equalDOF [expr int($nodeNum1-1+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+$Df_abut/$zSize4+($k-1)+($j-1)*$numZele4)] $nodeNum26 1 2 3 	     
            }}}}}}}

         if { $j == 1 && $k == [ expr int($numZeleabut+1)] } {
            set abutdeck_L $nodeNum26
            }
         set  nodeNum26  [expr $nodeNum26+1]   
  }
}

for { set j 1 } { $j <= [ expr $numYeleabut+1] } { incr j 1 } {
     equalDOF     [expr int($abutdeck_L+($j-1)*($numZeleabut+1))]     [expr int($DeckLnode+($j-1))]    1 2 3 4 5   
}
puts $meshFile "end coordinates"

for { set k 1 } { $k <= [ expr $numZeleabut] } { incr k 1 } {
	 for { set j 1 } { $j <= [ expr $numYeleabut] } { incr j 1 } {
	     set master [expr int($nodeNum25+($k-1))]
	     set slave  [expr int($master+$j*($numZeleabut+1))]
	     equalDOF    $master   $slave     1 2    
	     }
    }

puts $meshFile "Elements"

#####################################################
###### Develop Left Abutment Backwall Elements ######
#####################################################


set m 0
for {set i 1} {$i <= $numXeleabut} {incr i 1} {
  for {set j 1} {$j <= $numYeleabut} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
           
      set n1  [expr ($nodeNum25-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
      
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 44000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

puts $meshFile "end elements"


###################################################
###### Develop Right Abutment Backwall Nodes ######
###################################################

# Creates input file for GiD software
puts $meshFile "MESH shell dimension 3 ElemType Quadrilateral Nnode 4"
puts $meshFile "Coordinates"

set  nodeNum27  $nodeNum26

for { set j 1 } { $j <= [ expr $numYeleabut+1] } { incr j 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
    
         set xdim $L3

         if { $j <= [expr $numYele3 +1] } {
            set ydim [expr $Bzone1+$Bzone2+($j-1)*$ySize3]
            } else {
         if { $j == [expr $numYele3+1 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
            } else {
         if { $j == [expr $numYele3+2 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb/2.0] 
            } else {
         if { $j == [expr $numYele3+3 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb]                 
            } else {
         if { $j == [expr $numYele3+4 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
         if { $j == [expr $numYele3+5 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
             } else {
         if { $j == [expr $numYele3+6 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]                       
             } else {
         if { $j == [expr $numYele3+7 +1] } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]
             } else {
			 
         if { $j == [expr $numYele3+8+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+9+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+10+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]            
             } else {
         if { $j == [expr $numYele3+11+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+12+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+13+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]            
             } else {
         if { $j == [expr $numYele3+14+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]            
             } else {
         if { $j == [expr $numYele3+15+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile/2.0]            
             } else {
         if { $j == [expr $numYele3+16+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]            
             } else {			 
         if { $j == [expr $numYele3+17+1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]            
             } else {				 
         if { $j == [expr $numYele3+7+$numXele11+1+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]                                           
             } else {
         if { $j == [expr $numYele3+7+$numXele11+2+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]                                                                          
              } else {
         if { $j == [expr $numYele3+7+$numXele11+3+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]                                                                                                      
              } else {
         if { $j == [expr $numYele3+7+$numXele11+4+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]                                                                                                                                    
              } else {
         if { $j == [expr $numYele3+7+$numXele11+5+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]                                                                                                                                                                  
               } else {
         if { $j == [expr $numYele3+7+$numXele11+6+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
               } else {
         if { $j == [expr $numYele3+7+$numXele11+7+ 1]  } {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]            
               } else {
            set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb+($j-($numYele3+7+$numXele11+7+ 1))*$ySize5] 
               }}}}}}}}}}}}}}}}}}}}}}}}}

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum27 $xdim $ydim $zdim
         puts $meshFile "$nodeNum27 $xdim $ydim $zdim"
    
         fix $nodeNum27 0 0 0 0 0 1
         
         if { $k == 1 && $j == [expr int($numYele3+1+1)] } {
            equalDOF [expr int($np_g+$pile_elenum_abut)] $nodeNum27  1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+5+1)] } {
            equalDOF [expr int($np_g1+$pile_elenum_abut)] $nodeNum27 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+9+1)] } {
            equalDOF [expr int($np_g2+$pile_elenum_abut)] $nodeNum27 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+12+1)] } {
            equalDOF [expr int($np_g3+$pile_elenum_abut)] $nodeNum27 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+15+1)] } {
            equalDOF [expr int($np_g4+$pile_elenum_abut)] $nodeNum27 1 2 3 4 5 
            } else {
		 if { $k == 1 && $j == [expr int($numYele3+19+1)] } {
            equalDOF [expr int($np_g5+$pile_elenum_abut)] $nodeNum27 1 2 3 4 5 
            } else {
         if { $k == 1 && $j == [expr int($numYele3+7+$numXele11+6+1)] } {
            equalDOF [expr int($np_h+$pile_elenum_abut)] $nodeNum27  1 2 3 4 5 
            } else {
            equalDOF [expr int($nodeNum4-1+$Df_abut/$zSize4+($k-1)+($j-1)*$numZele4)] $nodeNum27   1 2 3
 	      }}}}}}}

         if { $j == 1 && $k == [ expr int($numZeleabut+1)] } {
            set abutdeck_R $nodeNum27
            }
         set  nodeNum27  [expr $nodeNum27+1]   
  }
}

for { set j 1 } { $j <= [ expr $numYeleabut+1] } { incr j 1 } {
     equalDOF    [expr int($abutdeck_R+($j-1)*($numZeleabut+1))]     [expr int($DeckRnode+($j-1))]     1 2 3 4 5
}
puts $meshFile "end coordinates"

for { set k 1 } { $k <= [ expr $numZeleabut] } { incr k 1 } {
	 for { set j 1 } { $j <= [ expr $numYeleabut] } { incr j 1 } {
	     set master [expr int($nodeNum26+($k-1))]
	     set slave  [expr int($master+$j*($numZeleabut+1))]
	     equalDOF    $master   $slave     1 2   
	     }
    }

puts $meshFile "Elements"

######################################################
###### Develop Right Abutment Backwall Elements ######
######################################################

set m 0
for {set i 1} {$i <= $numXeleabut} {incr i 1} {
  for {set j 1} {$j <= $numYeleabut} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
           
      set n1  [expr ($nodeNum26-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
      
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 44000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

puts $meshFile "end elements"

##################################################
###### Develop Left Abutment wingwall Nodes ######
##################################################

puts $meshFile "MESH shell dimension 3 ElemType Quadrilateral Nnode 4"
puts $meshFile "Coordinates"

set wingWele [expr $numXele3+2]

### Nodes for wing wall A #######

set  nodeNum28  $nodeNum27

for { set i 1 } { $i <= [ expr $wingWele] } { incr i 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
         
         if { $i <= [expr $numXele3+1] } {
            set xdim [expr $Lzone1+$Lzone2+($i-1)*$xSize3]
            } else {
            set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb]
            }

         set ydim [expr $Bzone1+$Bzone2] 

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum28 $xdim $ydim $zdim
         puts $meshFile "$nodeNum28 $xdim $ydim $zdim"
 
         fix $nodeNum28 0 0 0 0 0 1 
 
         equalDOF  [expr int($nodeNum1-1+($numXele_emb-$numXele3-2)*($numYele_deck+1)*$numZele4+$Df_abut/$zSize4+($k-1)+($i-1)*($numYele_deck+1)*$numZele4)]   $nodeNum28  1 2 3        

         set  nodeNum28  [expr $nodeNum28+1]   
  }
}

### Nodes for wing wall B #######

set  nodeNum29  $nodeNum28

for { set i 1 } { $i <= [ expr $wingWele] } { incr i 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
         
         if { $i <= [expr $numXele3+1] } {
            set xdim [expr $Lzone1+$Lzone2+($i-1)*$xSize3]
            } else {
            set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb]
            }

         set ydim [expr $B1+$Bzone5] 

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum29 $xdim $ydim $zdim
         puts $meshFile "$nodeNum29 $xdim $ydim $zdim"
 
         fix $nodeNum29 0 0 0 0 0 1
 
         equalDOF  [expr int($nodeNum1-1+($numXele_emb-$numXele3-2)*($numYele_deck+1)*$numZele4+($numYele_deck)*$numZele4+$Df_abut/$zSize4+($k-1)+($i-1)*($numYele_deck+1)*$numZele4)]   $nodeNum29  1 2 3        

         set  nodeNum29  [expr $nodeNum29+1]   
  }
}

###################################################
###### Develop Right Abutment wingwall Nodes ######
###################################################

### Nodes for wing wall C #######

set  nodeNum30  $nodeNum29

for { set i 1 } { $i <= [ expr $wingWele] } { incr i 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
         
         if { $i == 1 } {
            set xdim [expr $L3+$BP_emb/2.0]
            } else {
         if { $i == 2 } {
            set xdim [expr $L3+$BP_emb/2.0+$WInfemb]
            } else {
            set xdim [expr $L3+$BP_emb/2.0+$WInfemb+($i-2)*$xSize19]
            }}

         set ydim [expr $Bzone1+$Bzone2] 

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum30 $xdim $ydim $zdim
         puts $meshFile "$nodeNum30 $xdim $ydim $zdim"
 
         fix $nodeNum30 0 0 0 0 0 1
 
 
            equalDOF  [expr int($nodeNum4-1+($numYele_deck+1)*$numZele4+$Df_abut/$zSize4+($k-1)+($i-1)*($numYele_deck+1)*$numZele4)]   $nodeNum30  1 2 3        

 
         set  nodeNum30  [expr $nodeNum30+1]   
  }
}

### Nodes for wing wall D #######

set  nodeNum31  $nodeNum30

for { set i 1 } { $i <= [ expr $wingWele] } { incr i 1 } {
    for { set k 1 } { $k <= [ expr $numZeleabut+1] } { incr k 1 } {
         
         if { $i == 1 } {
            set xdim [expr $L3+$BP_emb/2.0]
            } else {
         if { $i == 2 } {
            set xdim [expr $L3+$BP_emb/2.0+$WInfemb]
            } else {
            set xdim [expr $L3+$BP_emb/2.0+$WInfemb+($i-2)*$xSize3]
            }}

         set ydim [expr $B1+$Bzone5] 

         set zdim [ expr $HZ+$Df_abut+($k-1)*$abutzsize]

         node  $nodeNum31 $xdim $ydim $zdim
         puts $meshFile "$nodeNum31 $xdim $ydim $zdim"
 
         fix $nodeNum31 0 0 0 0 0 1
 

         equalDOF  [expr int($nodeNum4-1+($numYele_deck+1)*$numZele4+($numYele_deck)*$numZele4+$Df_abut/$zSize4+($k-1)+($i-1)*($numYele_deck+1)*$numZele4)]   $nodeNum31  1 2 3        

 
         set  nodeNum31  [expr $nodeNum31+1]   
  }
}

puts $meshFile "end coordinates"
puts $meshFile "Elements"

#####################################################
###### Develop Left Abutment wingwall Elements ######
#####################################################

#### Elements for wing wall A #######
set m 0
for {set i 1} {$i <= $wingWele} {incr i 1} {
  for {set j 1} {$j <= 1} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
     if { $i != $wingWele } {      
      set n1  [expr ($nodeNum27-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    } else {
      set n1  [expr ($nodeNum27-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($nodeNum25+($k-1))] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    }  
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 54000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

#### Elements for wing wall B #######
set m 0
for {set i 1} {$i <= $wingWele} {incr i 1} {
  for {set j 1} {$j <= 1} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
     if { $i != $wingWele } {      
      set n1  [expr ($nodeNum28-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    } else {
      set n1  [expr ($nodeNum28-1)+$m*int($numZeleabut+1)+$k]
      set n2  [expr int($nodeNum25+$numYeleabut*($numZeleabut+1)+($k-1))] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    }  
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 54000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

######################################################
###### Develop Right Abutment wingwall Elements ######
######################################################
set m 0
for {set i 1} {$i <= $wingWele} {incr i 1} {
  for {set j 1} {$j <= 1} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
     if { $i != 1 } {      
      set n1  [expr ($nodeNum29-1)+($m-1)*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    } else {
      set n1  [expr int($nodeNum26+($k-1))]
      set n2  [expr ($nodeNum29-1)+$m*int($numZeleabut+1)+$k] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    }  
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 54000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

#### Elements for wing wall D #######
set m 0
for {set i 1} {$i <= $wingWele} {incr i 1} {
  for {set j 1} {$j <= 1} {incr j 1} {
    for {set k 1} {$k <= $numZeleabut} {incr k 1} {
     if { $i != 1 } {      
      set n1  [expr ($nodeNum30-1)+($m-1)*int($numZeleabut+1)+$k]
      set n2  [expr int($n1+$numZeleabut+1)] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    } else {
      set n1  [expr int($nodeNum26+$numYeleabut*($numZeleabut+1)+($k-1))]
      set n2  [expr ($nodeNum30-1)+$m*int($numZeleabut+1)+$k] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
    }  
     
     element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 54000
     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"

     set eleNum [expr $eleNum+1]
    }
    set m [expr ($m+1)]
  }
}

puts $meshFile "end elements"
