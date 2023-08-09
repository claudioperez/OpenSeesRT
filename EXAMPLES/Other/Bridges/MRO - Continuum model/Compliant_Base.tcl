
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates the nodes that are required for compliant base boundary.
# There is no user-defined parameter in the files but it is recommended to go through all the command lines before running it.
#========================================================================================


#################################################
#### Define Two Series of Nodes at Model Base ###
#################################################

model basic -ndm 3 -ndf 3

set ID   200000
set ID_2 400000

set nodeNum 0
for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {
    for {set k 1} { $k <= 1 } {incr k 1} {
	
	####### Coordinations in X direction
      if { $i <= [expr $numXele1+1] } {
	      set xdim [expr ($i-1)*$xSize1]
		 } else {
      if { $i > [expr $numXele1+1]  && $i <= [expr $numXele1+$numXele2+1] } { 		 
          set xdim [expr $Lzone1+($i-($numXele1+1))*$xSize2]
		 } else {
	  if { $i > [expr $numXele1+$numXele2+1]  && $i <= [expr $numXele1+$numXele2+$numXele3+1] } { 		 
          set xdim [expr $Lzone1+$Lzone2+($i-($numXele1+$numXele2+1))*$xSize3]
		 } else {	 	 
	  if { $i == [expr $numXele1+$numXele2+$numXele3+1+1] } {
          set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb]
		 } else {		 
	  if { $i > [expr $numXele1+$numXele2+$numXele3+1+1]  && $i <= [expr $numXele1+$numXele2+$numXele3+1+1+2] } {  
	     set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+($i-($numXele1+$numXele2+$numXele3+1+1))*$xSize4]
		 } else {		 
	  if { $i == [expr $numXele1+$numXele2+$numXele3+1+1+3] } {  
	     set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb+$WInfemb]  
		 } else {
	  if { $i > [expr $numXele1+$numXele2+$numXele3+4+1]   && $i <= [expr $numXele1+$numXele2+$numXele3+4+$numXele5+1]  } { 		 
          set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb+$WInfemb+($i-($numELE_a+1))*$xSize5]
		 } else {	  
	  if { $i > [expr $numELE_a+$numXele5+1]   && $i <= [expr $numELE_a+$numXele5+$numXele6+1]  } { 		 
          set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb+$WInfemb+$Lzone5+($i-($numELE_a+$numXele5+1))*$xSize6]
		 } else {  
	  if { $i > [expr $numELE_a+$numXele5+$numXele6+1]   && $i <= [expr $numELE_a+$numXele5+$numXele6+$numXele7+1]  } { 		 
          set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb+$WInfemb+$Lzone5+$Lzone6+($i-($numELE_a+$numXele5+$numXele6+1))*$xSize7]
		 } else { 
	  if { $i > [expr $numELE_b+1]   && $i <= [expr $numELE_b+$numXele8+1]  } { 		 
          set xdim [expr $L1+($i-($numELE_b+1))*$xSize8]
		 } else { 
	  if { $i > [expr $numELE_b+$numXele8+1]   && $i <= [expr $numELE_b+$numXele8+$numXele9+1]  } { 		 
          set xdim [expr $L1+$Lzone8+($i-($numELE_b+$numXele8+1))*$xSize9]
		 } else { 	
		 
	#### pile cap zone		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+1+1]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile]
		 } else { 	
	  if { $i > [expr $numELE_b+$numXele8+$numXele9+1+1]  && $i <= [expr $numELE_b+$numXele8+$numXele9+1+1+2]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+($i-($numELE_b+$numXele8+$numXele9+1+1))*$xSize10]
		 } else {	
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+5]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile]
		 } else { 	 		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+6]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11]
		 } else { 		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+7]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]
		 } else { 		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+8]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]
		 } else {		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+9]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]
		 } else {		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+10]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]
		 } else {				 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+11]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]
		 } else {			 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+12]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]
		 } else {				 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+13]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]
		 } else {		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+14]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]
		 } else {			 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+15]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]
		 } else {		 
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+1+1]  } { 	  
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]
		  } else {
	  if { $i > [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+1+1] && $i <= [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+3+1]   } { 		 
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+($i-($numELE_b+$numXele8+$numXele9+4+$numXele11+1+1))*$xSize10] 
		 } else {	  
	  if { $i == [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+1]  } { 	  
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	  
	     } else {
	  if { $i > [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+1] && $i <= [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+$numXele13+1] } { 	  
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+($i-($numELE_b+$numXele8+$numXele9+4+$numXele11+4+1))*$xSize13]	  
	     } else {	  
	  #########################
	  
	  if { $i > [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+$numXele13+1] && $i <= [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+$numXele13+$numXele14+1] } { 	  
          set xdim [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$Lzone13+($i-($numELE_b+$numXele8+$numXele9+4+$numXele11+4+$numXele13+1))*$xSize14]	  
	     } else {	
      if { $i > [expr $numELE_c+1] && $i <= [expr $numELE_c+$numXele15+1] } { 	  
          set xdim [expr $L2+($i-($numELE_c+1))*$xSize15]	  
	     } else {	   
      if { $i > [expr $numELE_c+$numXele15+1] && $i <= [expr $numELE_c+$numXele15+$numXele16+1] } { 	  
          set xdim [expr $L2+$Lzone15+($i-($numELE_c+$numXele15+1))*$xSize16]	  
	     } else {          
      if { $i > [expr $numELE_c+$numXele15+$numXele16+1] && $i <= [expr $numELE_c+$numXele15+$numXele16+$numXele17+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+($i-($numELE_c+$numXele15+$numXele16+1))*$xSize17]	  
	     } else {              
      if { $i == [expr $numELE_c+$numXele15+$numXele16+$numXele17+1+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb]	  
	     } else { 		 
      if { $i > [expr $numELE_c+$numXele15+$numXele16+$numXele17+1+1] && $i <= [expr $numELE_c+$numXele15+$numXele16+$numXele17+3+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+($i-($numELE_c+$numXele15+$numXele16+$numXele17+1+1))*$xSize18]	  
	     } else { 	
      if { $i == [expr $numELE_c+$numXele15+$numXele16+$numXele17+4+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+$BP_emb+$WInfemb]	  
	     } else { 		 
      if { $i > [expr $numELE_c+$numXele15+$numXele16+$numXele17+4+1] && $i <= [expr $numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+$BP_emb+$WInfemb+($i-($numELE_c+$numXele15+$numXele16+$numXele17+4+1))*$xSize19]	  
	     } else { 
      if { $i > [expr $numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+1] && $i <= [expr $numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+$numXele20+1] } { 	  
          set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+$BP_emb+$WInfemb+$Lzone19+($i-($numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+1))*$xSize20]	  
	     } else { 
      set xdim [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+$BP_emb+$WInfemb+$Lzone19+$Lzone20+($i-($numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+$numXele20+1))*$xSize21]
		}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}} 
		
	####### Coordinations in Y direction	
      if { $j <= [expr $numYele1+1] } {                                              
         set ydim [expr ($j-1)*$ySize1]       
         } else {
      if { $j > [expr $numYele1+1] && $j <= [expr ($numYele1+$numYele2+1) ] } {
         set ydim [expr $Bzone1+($j-($numYele1+1))*$ySize2]
         } else {		 
      if { $j > [expr ($numYele1+$numYele2+1)] && $j <= [expr ($numYele1+$numYele2+$numYele3+1) ] } {
         set ydim [expr $Bzone1+$Bzone2+($j-($numYele1+$numYele2+1))*$ySize3] 
         } else {
      if { $j == [expr ($numYele1+$numYele2+$numYele3+1+1) ] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
         } else {			 
	  if { $j > [expr ($numYele1+$numYele2+$numYele3+1+1)] && $j <= [expr ($numYele1+$numYele2+$numYele3+3+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+($j-($numYele1+$numYele2+$numYele3+1+1))*$ySize4]
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+4+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+5+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+6+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]	  
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]  	
         } else {
	############################	 
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+9)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+10)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]  	
         } else {		 
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+11)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+12)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+13)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]  		 
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+14)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {			 
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+15)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {			 
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+16)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]  	
         } else {	
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+17)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+18)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {	
	###################################	 
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+1+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]	
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+2+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]		 
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+3+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]		
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+4+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+5+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+6+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
         } else {
	  if { $j == [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+1)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
         } else {
	  if { $j > [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+1)] && $j <= [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+1)]} {
         set ydim [expr $B1+($j-($numYele1+$numYele2+$numYele3+7+$numXele11+7+1))*$ySize5]
         } else {
	  if { $j > [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+1)] && $j <= [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+$numYele6+1)]} {
         set ydim [expr $B1+$Bzone5+($j-($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+1))*$ySize6]		 
         } else {
	  if { $j > [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+$numYele6+1)] && $j <= [expr ($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+$numYele6+$numYele7+1)]} {
         set ydim [expr $B1+$Bzone5+$Bzone6+($j-($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+$numYele6+1))*$ySize7]
         }
		 }}}}}}}}}}}}}}}}}}}}}}}}}}}}
		 		 
    ####### Coordinations in Z direction
	  if { $k <= [expr ($numZele1+1)] } {
          set   zdim [expr ($k-1)*$zSize1]
         } else {
      if { $k > [expr ($numZele1+1)] && $k <= [expr ($numZele1+$numZele2+1)] } {
          set   zdim [expr $Czone1+($k-($numZele1+1))*$zSize2]
         } else {
      if { $k > [expr ($numZele1+$numZele2+1)] && $k <= [expr ($numZele1+$numZele2+$numZele3+1)] } {
          set   zdim [expr $Czone1+$Czone2+($k-($numZele1+$numZele2+1))*$zSize3]
         }}}
		
         set   nodeNum   [expr $nodeNum+1] 
		 
	   ## Top level
	   fix 	 [expr $ID+$m]    0 1 1 ; #User-defined --- release constraint in other directions if you want to apply motion in different directions.
	   ## Bottom level
	   fix 	 [expr $ID_2+$m]  1 1 1
	   ## Tie top level nodes to original base nodes
	   equalDOF [expr $ID+$m]  $nodeNum   1 2 3


   }
       set m [expr ($m+1)]
  } 
}


