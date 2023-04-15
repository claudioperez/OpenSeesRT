
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates the bridge deck per site condition at Meloland Road Overpass, in El Centro, CA.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================


################################
###### Develop Deck Nodes ######
################################

# Creates input file for GiD software
puts $meshFile "MESH beamcolumn dimension 3 ElemType Quadrilateral Nnode 4"
puts $meshFile "Coordinates"


set  nodeNum25  $nodeNum24
set x 1
for { set i 1 } { $i <= [ expr $numXele_D+1] } { incr i 1 } {
    for { set j 1 } { $j <= [ expr $numYele_D+1] } { incr j 1 } {
 
         if { $i <= [expr $Dzone1/$xSize_D1 +1] } {
            set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0+($i-1)*$xSize_D1]
            } else {
         if { $i > [expr $Dzone1/$xSize_D1+$Dzone2/$xSize_D2 +1] } {  
            set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0+$Ldeck-$Dzone3+($i-($Dzone1/$xSize_D1+$Dzone2/$xSize_D2 +1))*$xSize_D3]
            } else {
            set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0+$Dzone1+($i-($Dzone1/$xSize_D1 +1))*$xSize_D2]
            }}        
		  
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
         set zdim [expr $HZ+$Hpier]
         
         
         node  $nodeNum25 $xdim $ydim $zdim
         puts $meshFile "$nodeNum25 $xdim $ydim $zdim"
         
         if { $i == 1 && $j == 1 } {
         set DeckLnode $nodeNum25
            }
         if { $i == [expr $numXele_D+1] && $j == 1 } {
         set DeckRnode $nodeNum25
            }
            			
		 fix $nodeNum25 0 0 0 0 0 1
         
		 set  nodeNum25  [expr $nodeNum25+1]

    }
}

equalDOF     [expr int($nodeNum22+$Hpier/$psize+$Df/$zSize3)]   [expr int($nodeNum24-1+($numXele_D/2.0+1)*($numYele_D+1)-$numYele_D/2.0)]   1 2 3 4 5 6

for { set i 1 } { $i <= [ expr $numXele_D+1] } { incr i 1 } {
    for { set j 1 } { $j <=  $numYele_D } { incr j 1 } {
	     
        set master [expr int($nodeNum24+($i-1)*($numYele_D+1))]
        set slave  [expr $master+$j]
        equalDOF     $master    $slave  1 2  		 
              	   
       }
   }

puts $meshFile "end coordinates"
puts $meshFile "Elements"


###################################
###### Develop Deck Elements ######
###################################

set m 1
for {set i 1} {$i <= $numXele_D } {incr i 1} {
  for {set j 1} {$j <= $numYele_D } {incr j 1} {
      

      set n1  [expr int($nodeNum24+($i-1)*($numYele_D+1)+($j-1))]    
      set n2  [expr int($n1+($numYele_D+1))] 
      set n3  [expr int($n2+1)] 
      set n4  [expr int($n1+1)]
      
      element ShellMITC4 $eleNum $n1 $n2 $n3 $n4 32000

     puts $meshFile "$eleNum $n1 $n2 $n3 $n4 10"
     set eleNum [expr $eleNum+1] 
        
     }
   }  

puts $meshFile "end elements"
