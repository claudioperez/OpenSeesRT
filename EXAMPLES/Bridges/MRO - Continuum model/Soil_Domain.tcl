
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code creates soil domain per site condition at Meloland Road Overpass, in El Centro, CA.
# You need to make changes on the lines that the comment "User-Defined" appears. Please go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
#========================================================================================

#################################
###### Develop soil layers ######
#################################

# Creates input file for GiD software
puts $meshFile "MESH stdBrick dimension 3 ElemType Hexahedra Nnode 8"
puts $meshFile "Coordinates"

# Dimensions of piles and soil regions
set  BP_emb      0.32    					     ;#User-Defined --- pile diameter underneath the abutments
set  BP_pile     0.32   					     ;#User-Defined --- pile diameter underneath the pier
set  S_pier      3.32                            ;#User-Defined --- pile spacing underneath the pier
set  WInfpile    0.16							 ;#User-Defined --- soil element size at interface of piles under the pier
set  WInfemb     [expr $BP_emb/2.0]				 ;#User-Defined --- soil element size at interface of piles under the abutment
set  S_abut      [expr $S_pier+2.0*$WInfpile+2.0*$BP_pile] 		;#User-Defined --- pile spacing underneath the abutments
set  Df          2.0                                       		;#User-Defined --- Embedded depth of pile cap under the pier
set  Df_abut     3.0											;#User-Defined --- Embedded depth of pile cap under the abutment
set  LP_pier     15.50                                   	    ;#User-Defined --- pile length underneath the pier
set  LP_abut     18.0                                           ;#User-Defined --- pile length underneath the abutments
set  LP_abut     [expr $LP_abut-$Df_abut]

set  H_emb        6.0     			;#User-Defined --- Height of embankment
set  sloperatio   0.50    			;#User-Defined --- Side slope of embankment
set  sloperatio2  [expr 3.0/2.0]    ;#User-Defined --- Side slope of embankment in front of the abutment wall

# Defining zones 1 to 21 in the direction +X to optimize number of elements --- it starts from X=0 to X=L3
set   Lzone1      10.0		;#User-Defined 
set   Lzone2      4.0		;#User-Defined
set   Lzone3      3.0		;#User-Defined
set   Lzone5      [expr  ($H_emb-1*1.0)*$sloperatio2]
set   Lzone6      4.0		;#User-Defined
set   Lzone7      7.0		;#User-Defined
set   Lzone8      5.0		;#User-Defined
set   Lzone9      6.0		;#User-Defined
set   Lzone11     $S_pier
set   Lzone13     6.0		;#User-Defined
set   Lzone14     5.0		;#User-Defined
set   Lzone15     7.0		;#User-Defined
set   Lzone16     4.0		;#User-Defined
set   Lzone17     [expr  ($H_emb-1*1.0)*$sloperatio2]
set   Lzone19     3.0		;#User-Defined
set   Lzone20     4.0		;#User-Defined
set   Lzone21     10.0		;#User-Defined
set   L1   [expr $Lzone1+$Lzone2+$Lzone3+2.0*$WInfemb+$BP_emb+$Lzone5+$Lzone6+$Lzone7]
set   L2   [expr $L1+$Lzone8+$Lzone9+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$Lzone13+$Lzone14]
set   L3   [expr $L2+$Lzone15+$Lzone16+$Lzone17+$WInfemb+$BP_emb/2.0]

# Defining element sizes within zones 1 to 21 in the direction +X
set  xSize1    5.0		;#User-Defined --- Zone 1
set  xSize2    2.0		;#User-Defined --- Zone 2
set  xSize3    1.0		;#User-Defined --- Zone 3
set  xSize4    [expr $BP_emb/2.0]
set  xSize5    [expr $Lzone5/2.0]    
set  xSize6    $xSize2
set  xSize7    3.5		;#User-Defined --- Zone 7
set  xSize8    2.5		;#User-Defined --- Zone 8
set  xSize9    2.0		;#User-Defined --- Zone 9
set  xSize10   [expr $BP_pile/2.0]
set  xSize11   0.59		;#User-Defined --- Zone 11
set  xSize13   $xSize9
set  xSize14   $xSize8
set  xSize15   $xSize7
set  xSize16   2.0		;#User-Defined --- Zone 16
set  xSize17   [expr $Lzone17/2.0]
set  xSize18   [expr $BP_emb/2.0]	
set  xSize19   1.0		;#User-Defined --- Zone 19
set  xSize20   2.0		;#User-Defined --- Zone 20
set  xSize21   5.0		;#User-Defined --- Zone 21

# Number of elements within zones 1 to 21 in the direction +X
set  numXele1    [expr int($Lzone1/$xSize1)]
set  numXele2    [expr int($Lzone2/$xSize2)]
set  numXele3    [expr int($Lzone3/$xSize3)]
set  numELE_a    [expr $numXele1+$numXele2+$numXele3+4]
set  numXele5    [expr int($Lzone5/$xSize5)]
set  numXele6    [expr int($Lzone6/$xSize6)]
set  numXele7    [expr int($Lzone7/$xSize7)]
set  numXele8    [expr int($Lzone8/$xSize8)]
set  numXele9    [expr int($Lzone9/$xSize9)]
set  numXele11   10		;#User-Defined --- Zone 11 (within the pile cap under the brige pier)
set  numXele13   [expr int($Lzone13/$xSize13)]
set  numXele14   [expr int($Lzone14/$xSize14)]
set  numXele15   [expr int($Lzone15/$xSize15)]
set  numXele16   [expr int($Lzone16/$xSize16)]
set  numXele17   [expr int($Lzone17/$xSize17)]
set  numXele19   [expr int($Lzone19/$xSize19)]
set  numXele20   [expr int($Lzone20/$xSize20)]
set  numXele21   [expr int($Lzone21/$xSize21)]

set  numELE_b    [expr $numELE_a+$numXele5+$numXele6+$numXele7]
set  numELE_c    [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4+$numXele13+$numXele14]
set  numXELE_d     [expr $numELE_c+$numXele15+$numXele16+$numXele17+2]
set num_pile_X1  [expr $numXele1+$numXele2+$numXele3+1]
set num_pile_X2  [expr $numELE_b+$numXele8+$numXele9+1]
set num_pile_X3  [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+1]
set num_pile_X4  [expr $numXELE_d-1]


# Defining zones 1 to 7 in the direction +Y to optimize number of elements --- it starts from Y=0 to Y=B1
set Bzone1 10.0			;#User-Defined --- Zone 1
set Bzone2 [expr $H_emb/$sloperatio]
set Bzone3 [expr (10.0-2.0*$WInfemb-2.0*$BP_emb-2.0*$WInfpile-$S_abut)/2.0] ; # User-Defined --- difference between width of the embankment width and pile cap width under the bridge pier
set Bzone5 $Bzone3
set Bzone6 [expr $H_emb/$sloperatio]
set Bzone7 10.0			;#User-Defined --- Zone 7

# Defining element sizes within zones 1 to 7 in the direction +Y
set  ySize1   5.0		;#User-Defined
set  ySize2   [expr $Bzone2/2.0]
set  ySize3   $Bzone3
set  ySize4   [expr $BP_emb/2.0]
set  ySize5   $Bzone5
set  ySize6   [expr $Bzone6/2.0]
set  ySize7   5.0		;#User-Defined

# Number of elements within zones 1 to 7 in the direction +Y
set numYele1  [expr int($Bzone1/$ySize1)]
set numYele2  [expr int($Bzone2/$ySize2)]
set numYele3  [expr int($Bzone3/$ySize3)]
set numYele5  [expr int($Bzone5/$ySize5)]
set numYele6  [expr int($Bzone6/$ySize6)]
set numYele7  [expr int($Bzone7/$ySize7)]

# Parameters to be used for locating specific soil nodes and elements
set B1 [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
set numYele_dum [expr $numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5]
set num_pile_Y1  [expr $numYele1+$numYele2+$numYele3+1]
set num_pile_Y2  [expr $numYele1+$numYele2+$numYele3+4]
set num_pile_Y3  [expr $numYele1+$numYele2+$numYele3+7+$numXele11+1]
set num_pile_Y4  [expr $numYele1+$numYele2+$numYele3+7+$numXele11+4]


# Defining zones 1 to 4 in the direction +Z --- it starts from Z=0 to Z=HZ
set Czone1 5.0 		;#User-Defined
set Czone2 8.0 		;#User-Defined
set Czone3 7.0 		;#User-Defined
set Czone4 $H_emb   ;#User-Defined --- Height of embankment
set HZ      [expr $Czone1+$Czone2+$Czone3]

# Element size in direction +Z
set zSize1  2.5
set zSize2  2.0
set zSize3  1.0
set zSize4  1.0

# Number of soil elements in direction +Z
set numZele1  [expr int($Czone1/$zSize1)]
set numZele2  [expr int($Czone2/$zSize2)]
set numZele3  [expr int($Czone3/$zSize3)]
set numZele4  [expr int($Czone4/$zSize4)]

# Total number of elements
set numXele    [expr int($numELE_c+$numXele15+$numXele16+$numXele17+4+$numXele19+$numXele20+$numXele21)]
set numYele    [expr int($numYele1+$numYele2+$numYele3+7+$numXele11+7+$numYele5+$numYele6+$numYele7)]
set numZele    [expr int($numZele1+$numZele2+$numZele3)]
set   numYele_deck [expr $Bzone3/$ySize3+7+$numXele11+7+$Bzone5/$ySize5]

# Total number of nodes
set   numXnode   [expr $numXele+1]
set   numYnode   [expr $numYele+1]
set   numZnode   [expr $numZele+1]

# To be used to locate specific elements and nodes
set   numZele_a 	[expr $numZele1+$numZele2+$numZele3+$numZele4]
set   numZele_b  	[expr $numZele1+$numZele2+$numZele3]
set   numXele_emb   [expr $numXele1+$numXele2+$numXele3+2]
set   numXnode_emb  [expr $numXele_emb+1]
set   num_pile_Z1   [expr int(($HZ-$LP_abut)/$zSize1)]
set   num_pile_Z2   [expr int(($HZ-($LP_pier+$Df))/$zSize1)]


####################################
###### Develop soil Materials ######
####################################

##layer 1
set Bs1      300000.0; #User-Defined
set Gs1      60000.0; #User-Defined
set Es1     [expr 9.0*$Bs1*$Gs1/(3.0*$Bs1+$Gs1)]
set nou1    0.3; #User-Defined
set rhos1   1.5; #User-Defined
set cohes1  35.9; #User-Defined
set peakShearStra 0.1; #User-Defined
nDMaterial PressureIndependMultiYield 1 3 $rhos1  $Gs1 $Bs1 $cohes1 $peakShearStra

##layer 2
set B2            200000.0; #User-Defined
set G2            75000.0; #User-Defined
set Es2           [expr 9.0*$B2*$G2/(3.0*$B2+$G2)]
set nou2           0.3; #User-Defined
set rhos2          1.9; #User-Defined
set Phi2           33.0 ; #User-Defined
set cohes2         76.6 ; #User-Defined
set refPress       80.0; #User-Defined
set pressDependCoe 0.5; #User-Defined
set PTAng          27.0; #User-Defined
set contrac        0.07; #User-Defined
set dilat1         0.40; #User-Defined
set dilat2         2.0; #User-Defined
set liquefac1      0.0; #User-Defined
set liquefac2      0.0; #User-Defined
set liquefac3      0.0; #User-Defined
nDMaterial PressureDependMultiYield 2 3 $rhos2  $G2 $B2 $Phi2 $peakShearStra $refPress $pressDependCoe $PTAng $contrac $dilat1 $dilat2 $liquefac1 $liquefac2 $liquefac3 

##layer 3
set B3      750000.0; #User-Defined
set G3      150000.0; #User-Defined
set Es3     [expr 9.0*$B3*$G3/(3.0*$B3+$G3)]
set nou3    0.3; #User-Defined
set rhos3   1.8; #User-Defined
set cohes3  76.6; #User-Defined
nDMaterial PressureIndependMultiYield 3 3 $rhos3  $G3 $B3 $cohes3 $peakShearStra

##layer 4
set B4      200000.0; #User-Defined
set G4      75000.0; #User-Defined
set Es4     [expr 9.0*$B4*$G4/(3.0*$B4+$G4)]
set nou4    0.3; #User-Defined
set rhos4   1.9; #User-Defined
set cohes4   86.2; #User-Defined
set Phi4           33.0 ; #User-Defined 
set refpress       80.0; #User-Defined
set pressDependCoe 0.5; #User-Defined
set PTAng          27.0; #User-Defined
set contract       0.05; #User-Defined
set dilat1         0.60; #User-Defined
set dilat2         3.0; #User-Defined
set liquefac1      0.0; #User-Defined
set liquefac2      0.0; #User-Defined
set liquefac3      0.0; #User-Defined
nDMaterial PressureDependMultiYield 4 3 $rhos4  $G4 $B4 $Phi4 $peakShearStra $refPress $pressDependCoe $PTAng $contrac $dilat1 $dilat2 $liquefac1 $liquefac2 $liquefac3 

##layer 5
set B5      750000.0; #User-Defined
set G5      150000.0; #User-Defined
set Es5     [expr 9.0*$B5*$G5/(3.0*$B5+$G5)]
set nou5    0.3; #User-Defined
set rhos5   1.8; #User-Defined
set cohes5  86.2; #User-Defined
nDMaterial PressureIndependMultiYield 5 3 $rhos5  $G5 $B5 $cohes5 $peakShearStra

##Pile Cap
set Ecap       22.0e+6; #User-Defined
set noucap     0.20;    #User-Defined
set rhocap     2.40;    #User-Defined
set Gcap       [expr $Ecap/2.0/(1+0.2)]
set Bcap       [expr $Ecap/3.0/(1-2*0.2)]
set peakShearStracap  1.0; #User-Defined
set cohescap 200.0;        #User-Defined
nDMaterial ElasticIsotropic 6 $Ecap $noucap $rhocap

## Embankment
set nouemb      0.4;     #User-Defined
set Gemb        19400.0; #User-Defined
set Bemb        90000.0; #User-Defined
set Esemb       [expr $Gemb*2.0*(1+$nouemb)]
set rhosemb     1.6;     #User-Defined
set cohesemb    20.0;    #User-Defined
nDMaterial PressureIndependMultiYield 7 3 $rhosemb  $Gemb $Bemb $cohesemb $peakShearStra

###################################
###### Define gravity forces ######
###################################
set gravX  0.0000
set gravY  0.0000
set gravZ -9.8100


##########################
###### Create Nodes ######
##########################

set nodeNum 0
for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {
    for {set k 1} { $k <= [expr $numZnode] } {incr k 1} {
	
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
         node  $nodeNum  $xdim $ydim $zdim
         puts $meshFile "$nodeNum  $xdim $ydim $zdim"

	 ############## To locate pile nodes under left abutment wall	 
     ### PILE A ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y1+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a $nodeNum
        set p2_a [expr $p1_a+$numYnode*$numZnode] 
        set p3_a [expr $p1_a+2*$numYnode*$numZnode]
        }
	 ### PILE A1 ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a1 $nodeNum
        set p2_a1 [expr $p1_a1+$numYnode*$numZnode] 
        set p3_a1 [expr $p1_a1+2*$numYnode*$numZnode]
        }
	 ### PILE A2 ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a2 $nodeNum
        set p2_a2 [expr $p1_a2+$numYnode*$numZnode] 
        set p3_a2 [expr $p1_a2+2*$numYnode*$numZnode]
        }
	 ### PILE A3 ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a3 $nodeNum
        set p2_a3 [expr $p1_a3+$numYnode*$numZnode] 
        set p3_a3 [expr $p1_a3+2*$numYnode*$numZnode]
        }
	 ### PILE A4 ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a4 $nodeNum
        set p2_a4 [expr $p1_a4+$numYnode*$numZnode] 
        set p3_a4 [expr $p1_a4+2*$numYnode*$numZnode]
        }
	 ### PILE A5 ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z1+1] } { 
        set p1_a5 $nodeNum
        set p2_a5 [expr $p1_a5+$numYnode*$numZnode] 
        set p3_a5 [expr $p1_a5+2*$numYnode*$numZnode]
        }				
     ### PILE B ###    
     if { $i == [expr $num_pile_X1+1] && $j == [expr $num_pile_Y4+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_b $nodeNum
        set p2_b [expr $p1_b+$numYnode*$numZnode] 
        set p3_b [expr $p1_b+2*$numYnode*$numZnode]
        }
		
	 ############## To locate pile nodes under left abutment wall			
     ### PILE 1 ###    
     if { $i == [expr $num_pile_X2+1] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z2+1] } { 
        set p1_1 $nodeNum
        set p2_1 [expr $p1_1+$numYnode*$numZnode] 
        set p3_1 [expr $p1_1+2*$numYnode*$numZnode]
        }
	 ### PILE 2 ###    
     if { $i == [expr $num_pile_X2+1] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z2+1] } { 
        set p1_2 $nodeNum
        set p2_2 [expr $p1_2+$numYnode*$numZnode] 
        set p3_2 [expr $p1_2+2*$numYnode*$numZnode]
        }
	 ### PILE 3 ###    
     if { $i == [expr $num_pile_X2+1] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z2+1] } { 
        set p1_3 $nodeNum
        set p2_3 [expr $p1_3+$numYnode*$numZnode] 
        set p3_3 [expr $p1_3+2*$numYnode*$numZnode]
        }		
	 ### PILE 4 ###    
     if { $i == [expr $num_pile_X2+1] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z2+1] } { 
        set p1_4 $nodeNum
        set p2_4 [expr $p1_4+$numYnode*$numZnode] 
        set p3_4 [expr $p1_4+2*$numYnode*$numZnode]
        }					
     ### PILE 5 ###    
     if { $i == [expr $num_pile_X2+1] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z2+1] } { 
        set p1_5 $nodeNum
        set p2_5 [expr $p1_5+$numYnode*$numZnode] 
        set p3_5 [expr $p1_5+2*$numYnode*$numZnode]
        }				
     ### PILE 6 ###    
     if { $i == [expr $num_pile_X2+5] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z2+1] } { 
        set p1_6 $nodeNum
        set p2_6 [expr $p1_6+$numYnode*$numZnode] 
        set p3_6 [expr $p1_6+2*$numYnode*$numZnode]
        }
	 ### PILE 7 ###    
     if { $i == [expr $num_pile_X2+5] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z2+1] } { 
        set p1_7 $nodeNum
        set p2_7 [expr $p1_7+$numYnode*$numZnode] 
        set p3_7 [expr $p1_7+2*$numYnode*$numZnode]
        }
	 ### PILE 8 ###    
     if { $i == [expr $num_pile_X2+5] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z2+1] } { 
        set p1_8 $nodeNum
        set p2_8 [expr $p1_8+$numYnode*$numZnode] 
        set p3_8 [expr $p1_8+2*$numYnode*$numZnode]
        }		
	 ### PILE 9 ###    
     if { $i == [expr $num_pile_X2+5] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z2+1] } { 
        set p1_9 $nodeNum
        set p2_9 [expr $p1_9+$numYnode*$numZnode] 
        set p3_9 [expr $p1_9+2*$numYnode*$numZnode]
        }					
     ### PILE 10 ###    
     if { $i == [expr $num_pile_X2+5] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z2+1] } { 
        set p1_10 $nodeNum
        set p2_10 [expr $p1_10+$numYnode*$numZnode] 
        set p3_10 [expr $p1_10+2*$numYnode*$numZnode]
        }
    ### PILE 11 ###    
     if { $i == [expr $num_pile_X2+8] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z2+1] } { 
        set p1_11 $nodeNum
        set p2_11 [expr $p1_11+$numYnode*$numZnode] 
        set p3_11 [expr $p1_11+2*$numYnode*$numZnode]
        }
	 ### PILE 12 ###    
     if { $i == [expr $num_pile_X2+8] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z2+1] } { 
        set p1_12 $nodeNum
        set p2_12 [expr $p1_12+$numYnode*$numZnode] 
        set p3_12 [expr $p1_12+2*$numYnode*$numZnode]
        }
	 ### PILE 13 ###    
     if { $i == [expr $num_pile_X2+8] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z2+1] } { 
        set p1_13 $nodeNum
        set p2_13 [expr $p1_13+$numYnode*$numZnode] 
        set p3_13 [expr $p1_13+2*$numYnode*$numZnode]
        }		
	 ### PILE 14 ###    
     if { $i == [expr $num_pile_X2+8] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z2+1] } { 
        set p1_14 $nodeNum
        set p2_14 [expr $p1_14+$numYnode*$numZnode] 
        set p3_14 [expr $p1_14+2*$numYnode*$numZnode]
        }					
     ### PILE 15 ###    
     if { $i == [expr $num_pile_X2+8] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z2+1] } { 
        set p1_15 $nodeNum
        set p2_15 [expr $p1_15+$numYnode*$numZnode] 
        set p3_15 [expr $p1_15+2*$numYnode*$numZnode]
        }

	 ### PILE 16 ###    
     if { $i == [expr $num_pile_X2+11] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z2+1] } { 
        set p1_16 $nodeNum
        set p2_16 [expr $p1_16+$numYnode*$numZnode] 
        set p3_16 [expr $p1_16+2*$numYnode*$numZnode]
        }
	 ### PILE 17 ###    
     if { $i == [expr $num_pile_X2+11] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z2+1] } { 
        set p1_17 $nodeNum
        set p2_17 [expr $p1_17+$numYnode*$numZnode] 
        set p3_17 [expr $p1_17+2*$numYnode*$numZnode]
        }
	 ### PILE 18 ###    
     if { $i == [expr $num_pile_X2+11] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z2+1] } { 
        set p1_18 $nodeNum
        set p2_18 [expr $p1_18+$numYnode*$numZnode] 
        set p3_18 [expr $p1_18+2*$numYnode*$numZnode]
        }		
	 ### PILE 19 ###    
     if { $i == [expr $num_pile_X2+11] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z2+1] } { 
        set p1_19 $nodeNum
        set p2_19 [expr $p1_19+$numYnode*$numZnode] 
        set p3_19 [expr $p1_19+2*$numYnode*$numZnode]
        }					
     ### PILE 20 ###    
     if { $i == [expr $num_pile_X2+11] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z2+1] } { 
        set p1_20 $nodeNum
        set p2_20 [expr $p1_20+$numYnode*$numZnode] 
        set p3_20 [expr $p1_20+2*$numYnode*$numZnode]
        }

	 ### PILE 21 ###    
     if { $i == [expr $num_pile_X2+15] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z2+1] } { 
        set p1_21 $nodeNum
        set p2_21 [expr $p1_21+$numYnode*$numZnode] 
        set p3_21 [expr $p1_21+2*$numYnode*$numZnode]
        }
	 ### PILE 22 ###    
     if { $i == [expr $num_pile_X2+15] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z2+1] } { 
        set p1_22 $nodeNum
        set p2_22 [expr $p1_22+$numYnode*$numZnode] 
        set p3_22 [expr $p1_22+2*$numYnode*$numZnode]
        }
	 ### PILE 23 ###    
     if { $i == [expr $num_pile_X2+15] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z2+1] } { 
        set p1_23 $nodeNum
        set p2_23 [expr $p1_23+$numYnode*$numZnode] 
        set p3_23 [expr $p1_23+2*$numYnode*$numZnode]
        }		
	 ### PILE 24 ###    
     if { $i == [expr $num_pile_X2+15] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z2+1] } { 
        set p1_24 $nodeNum
        set p2_24 [expr $p1_24+$numYnode*$numZnode] 
        set p3_24 [expr $p1_24+2*$numYnode*$numZnode]
        }					
     ### PILE 25 ###    
     if { $i == [expr $num_pile_X2+15] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z2+1] } { 
        set p1_25 $nodeNum
        set p2_25 [expr $p1_25+$numYnode*$numZnode] 
        set p3_25 [expr $p1_25+2*$numYnode*$numZnode]
        }
		
	 ############## To locate pile nodes under right abutment wall			
     ### PILE G ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y1+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g $nodeNum
        set p2_g [expr $p1_g+$numYnode*$numZnode] 
        set p3_g [expr $p1_g+2*$numYnode*$numZnode]
        }
		
	 ### PILE G1 ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y2+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g1 $nodeNum
        set p2_g1 [expr $p1_g1+$numYnode*$numZnode] 
        set p3_g1 [expr $p1_g1+2*$numYnode*$numZnode]
        }
	 ### PILE G2 ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y2+5] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g2 $nodeNum
        set p2_g2 [expr $p1_g2+$numYnode*$numZnode] 
        set p3_g2 [expr $p1_g2+2*$numYnode*$numZnode]
        }
	 ### PILE G3 ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y2+8] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g3 $nodeNum
        set p2_g3 [expr $p1_g3+$numYnode*$numZnode] 
        set p3_g3 [expr $p1_g3+2*$numYnode*$numZnode]
        }
	 ### PILE G4 ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y2+11] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g4 $nodeNum
        set p2_g4 [expr $p1_g4+$numYnode*$numZnode] 
        set p3_g4 [expr $p1_g4+2*$numYnode*$numZnode]
        }
	 ### PILE G5 ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y2+15] && $k == [expr $num_pile_Z1+1] } { 
        set p1_g5 $nodeNum
        set p2_g5 [expr $p1_g5+$numYnode*$numZnode] 
        set p3_g5 [expr $p1_g5+2*$numYnode*$numZnode]
        }		
		
     ### PILE H ###    
     if { $i == [expr $num_pile_X4+1] && $j == [expr $num_pile_Y4+1] && $k == [expr $num_pile_Z1+1] } { 
        set p1_h $nodeNum
        set p2_h [expr $p1_h+$numYnode*$numZnode] 
        set p3_h [expr $p1_h+2*$numYnode*$numZnode]
        }

     ### Pier ###
     if { $i == [expr int($num_pile_X2+3+$numXele11/2.0+1)] && $j == [expr int($num_pile_Y1+6+$numXele11/2.0+1)] && $k == [expr int($numZnode-1-$Df/$zSize3+1)] } { 
        set pierbase $nodeNum
        }

   }
  } 
}


##############################################
###### Define soil nodes at embankments ######
##############################################

# Embankment on the left side
set   nodeNum1   [expr $nodeNum+1]

for {set i 1} {$i <= $numXnode_emb} {incr i 1} {
set m 1
  for {set j 1} {$j <= $numZele4} {incr j 1} {
  
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
	      if { $i > [expr $numXele1+$numXele2+$numXele3+1+1]  } {  
	         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0]
		     } 
             }}}}
                 
          if { $j != $numZele4 } {
             set ydim [expr $Bzone1+($zSize4/$sloperatio)*$j]
             set ydim_2 [expr $ydim+($Bzone2-($zSize4/$sloperatio)*$j)/2.0]
             } else {
             set ydim [expr $Bzone1+($zSize4/$sloperatio)*$j-($zSize4/$sloperatio)/2.0]
             set ydim_2 [expr $ydim+($Bzone2-($zSize4/$sloperatio)*$j)/2.0+($zSize4/$sloperatio)/4.0]
             }
                         
             set zdim   [expr $Czone1+$Czone2+$Czone3+$m*$zSize4]
             set m [expr $m+1]
             
             node  $nodeNum1  $xdim $ydim $zdim
             puts $meshFile "$nodeNum1  $xdim $ydim $zdim"
             set   nodeNum1   [expr $nodeNum1+1]
             
             node  $nodeNum1  $xdim $ydim_2 $zdim 
             puts $meshFile "$nodeNum1  $xdim $ydim_2 $zdim"
             set   nodeNum1   [expr $nodeNum1+1]
             
      }
   }    

set   nodeNum2   $nodeNum1

for { set i 1 } { $i <= $numXnode_emb } { incr i 1 } {
  for {set j 1} {$j <=  [expr $numYele_deck+1] } {incr j 1} {
      for {set k 1} { $k <= $numZele4 } {incr k 1} {
      
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
	      if { $i >  [expr $numXele1+$numXele2+$numXele3+1+1] } {  
	         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0]
		     } 
             }}}}
                        
          if { $j <= [expr $numYele3+1] } {
             set ydim [expr $Bzone1+$Bzone2+(($j-1)*$ySize3)]
             } else {
          if { $j == [expr $numYele3+1+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
             } else {			 
	      if { $j > [expr $numYele3+1+1] && $j <= [expr $numYele3+1+1+2] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+($j-($numYele3+1+1))*$ySize4]
             } else {
	      if { $j == [expr ($numYele3+1+1+3)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
	      if { $j == [expr ($numYele3+1+1+4)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
            } else {
	      if { $j == [expr ($numYele3+1+1+5)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]	  
             } else {
	      if { $j == [expr ($numYele3+1+1+6)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]  	
             } else {
			 
      ############################	 
	  if { $j == [expr ($numYele3+9)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele3+10)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]  	
         } else {		 
	  if { $j == [expr ($numYele3+11)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele3+12)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele3+13)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]  		 
         } else {
	  if { $j == [expr ($numYele3+14)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {			 
	  if { $j == [expr ($numYele3+15)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {			 
	  if { $j == [expr ($numYele3+16)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]  	
         } else {	
	  if { $j == [expr ($numYele3+17)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele3+18)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {	
	  ###################################	 
	
	      if { $j == [expr ($numYele3+7+$numXele11+1+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]	
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+2+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]		 
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+3+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]		
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+4+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	
             } else {
          if { $j == [expr ($numYele3+7+$numXele11+5+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+6+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
             } else {
	      if { $j == [expr $numYele3+7+$numXele11+7+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
             } else {
             set ydim [expr $B1+($j-($numYele3+7+$numXele11+7+1))*$ySize5]
             }}}}}}}}}}}}}}}}}}}}}}}}            
                  
             set zdim   [expr $Czone1+$Czone2+$Czone3+$k*$zSize4]
      
             node  $nodeNum2  $xdim $ydim $zdim
             puts $meshFile "$nodeNum2  $xdim $ydim $zdim"
             set   nodeNum2   [expr $nodeNum2+1]
             
         }
      }
   }
   
set   nodeNum3    $nodeNum2

for {set i 1} {$i <= $numXnode_emb} {incr i 1} {
  set m 1
  for {set j 1} {$j <= $numZele4 } {incr j 1} {
      
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
	      if { $i >  [expr $numXele1+$numXele2+$numXele3+1+1] } {  
	         set xdim [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0]
		     } 
             }}}} 
                         
          if { $j != 1 } {
             set ydim [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*($numZele4-$j+1)]
             } else {
             set ydim [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*($numZele4-$j+1)+($zSize4/$sloperatio)/2.0]
             }
             set ydim_2 [expr $ydim-($ydim-($B1+$Bzone5))/2.0]
             
             set zdim   [expr $Czone1+$Czone2+$Czone3+($Czone4-($m-1)*$zSize4)]
             set m [expr $m+1]     
             
		  
             node  $nodeNum3  $xdim $ydim $zdim
             puts $meshFile "$nodeNum3  $xdim $ydim $zdim"
             set   nodeNum3   [expr $nodeNum3+1]
             
             node  $nodeNum3  $xdim $ydim_2 $zdim 
             puts $meshFile "$nodeNum3  $xdim $ydim_2 $zdim"
             set   nodeNum3   [expr $nodeNum3+1]
            
      }
   }

# Embankment on the right side

set   nodeNum4  $nodeNum3 

for {set i 1} {$i <= $numXnode_emb} {incr i 1} {
set m 1
  for {set j 1} {$j <= $numZele4} {incr j 1} {
    
          if { $i == 1 } {
	         set xdim [expr $L3]
		     } else {
          if { $i  == 2 } { 		 
             set xdim [expr $L3+$BP_emb/2.0]
		     } else {
          if { $i == 3 } { 		 
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb]
		     } else {	 	 
          if { $i > 3 && $i <= [expr 2+$numXele19+1] } {
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb+($i-(3))*$xSize19]
		     } else {		 
          if { $i > [expr 2+$numXele19+1]  &&  $i <= [expr 2+$numXele19+$numXele20+1]} {  
	         set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+($i-(2+$numXele19+1))*$xSize20]
		     } else {
			 set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+$Lzone20+($i-(2+$numXele19+$numXele20+1))*$xSize21]
             }}}}}
                  
          if { $j != $numZele4 } {
             set ydim [expr $Bzone1+($zSize4/$sloperatio)*$j]
             set ydim_2 [expr $ydim+($Bzone2-($zSize4/$sloperatio)*$j)/2.0]
             } else {
             set ydim [expr $Bzone1+($zSize4/$sloperatio)*$j-($zSize4/$sloperatio)/2.0]
             set ydim_2 [expr $ydim+($Bzone2-($zSize4/$sloperatio)*$j)/2.0+($zSize4/$sloperatio)/4.0]
             }
                         
             set zdim   [expr $Czone1+$Czone2+$Czone3+$m*$zSize4]
             set m [expr $m+1]            
		 
             node  $nodeNum4  $xdim $ydim $zdim
             puts $meshFile "$nodeNum4  $xdim $ydim $zdim"
             set   nodeNum4   [expr $nodeNum4+1]
             
             node  $nodeNum4  $xdim $ydim_2 $zdim 
             puts $meshFile "$nodeNum4  $xdim $ydim_2 $zdim"
             set   nodeNum4   [expr $nodeNum4+1]
        
      }
   }    

set   nodeNum5   $nodeNum4

for { set i 1 } { $i <= $numXnode_emb } { incr i 1 } {
  for {set j 1} {$j <=  [expr $numYele_deck+1] } {incr j 1} {
      for {set k 1} { $k <= $numZele4 } {incr k 1} {
      
          if { $i == 1 } {
	         set xdim [expr $L3]
		     } else {
          if { $i  == 2 } { 		 
             set xdim [expr $L3+$BP_emb/2.0]
		     } else {
	      if { $i == 3 } { 		 
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb]
		     } else {	 	 
	      if { $i > 3 && $i <= [expr 2+$numXele19+1] } {
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb+($i-(3))*$xSize19]
		     } else {		 
	      if { $i > [expr 2+$numXele19+1]  &&  $i <= [expr 2+$numXele19+$numXele20+1]} {  
	         set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+($i-(2+$numXele19+1))*$xSize20]
		     } else {
			 set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+$Lzone20+($i-(2+$numXele19+$numXele20+1))*$xSize21]
             }}}}}
                        
          if { $j <= [expr $numYele3+1] } {
             set ydim [expr $Bzone1+$Bzone2+(($j-1)*$ySize3)]
             } else {
          if { $j == [expr $numYele3+1+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
             } else {			 
	      if { $j > [expr $numYele3+1+1] && $j <= [expr $numYele3+1+1+2] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+($j-($numYele3+1+1))*$ySize4]
             } else {
	      if { $j == [expr ($numYele3+1+1+3)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
	      if { $j == [expr ($numYele3+1+1+4)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
            } else {
	      if { $j == [expr ($numYele3+1+1+5)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]	  
             } else {
	      if { $j == [expr ($numYele3+1+1+6)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]  	
             } else {
			 
	  ############################	 
	  if { $j == [expr ($numYele3+9)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele3+10)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]  	
         } else {		 
	  if { $j == [expr ($numYele3+11)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele3+12)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]  	
         } else {
	  if { $j == [expr ($numYele3+13)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]  		 
         } else {
	  if { $j == [expr ($numYele3+14)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {			 
	  if { $j == [expr ($numYele3+15)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {			 
	  if { $j == [expr ($numYele3+16)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]  	
         } else {	
	  if { $j == [expr ($numYele3+17)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {	
	  if { $j == [expr ($numYele3+18)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {	
	###################################	 
	
	      if { $j == [expr ($numYele3+7+$numXele11+1+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]	
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+2+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]		 
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+3+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]		
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+4+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	
             } else {
          if { $j == [expr ($numYele3+7+$numXele11+5+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]
             } else {
	      if { $j == [expr ($numYele3+7+$numXele11+6+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
             } else {
	      if { $j == [expr $numYele3+7+$numXele11+7+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
             } else {
             set ydim [expr $B1+($j-($numYele3+7+$numXele11+7+1))*$ySize5]
             }}}}}}}}}}}}}}}}}}}}}}}}            
                        
             set zdim   [expr $Czone1+$Czone2+$Czone3+$k*$zSize4]
  
             node  $nodeNum5  $xdim $ydim $zdim
             puts $meshFile "$nodeNum5  $xdim $ydim $zdim"
             set   nodeNum5  [expr $nodeNum5+1]
             
         }
      }
   } 
set   nodeNum6    $nodeNum5

for {set i 1} {$i <= $numXnode_emb} {incr i 1} {
  set m 1
  for {set j 1} {$j <= $numZele4 } {incr j 1} {
      
          if { $i == 1 } {
	         set xdim [expr $L3]
		     } else {
          if { $i  == 2 } { 		 
             set xdim [expr $L3+$BP_emb/2.0]
		     } else {
	      if { $i == 3 } { 		 
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb]
		     } else {	 	 
	      if { $i > 3 && $i <= [expr 2+$numXele19+1] } {
             set xdim [expr $L3+$BP_emb/2.0+$WInfemb+($i-(3))*$xSize19]
		     } else {		 
	      if { $i > [expr 2+$numXele19+1]  &&  $i <= [expr 2+$numXele19+$numXele20+1]} {  
	         set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+($i-(2+$numXele19+1))*$xSize20]
		     } else {
			 set xdim [expr $L3+$BP_emb/2.0+$WInfemb+$Lzone19+$Lzone20+($i-(2+$numXele19+$numXele20+1))*$xSize21]
             }}}}}            
             
          if { $j != 1 } {
             set ydim [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*($numZele4-$j+1)]
             } else {
             set ydim [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*($numZele4-$j+1)+($zSize4/$sloperatio)/2.0]
             }
             set ydim_2 [expr $ydim-($ydim-($B1+$Bzone5))/2.0]
             
             set zdim   [expr $Czone1+$Czone2+$Czone3+($Czone4-($m-1)*$zSize4)]
             set m [expr $m+1]     
             
		     
             node  $nodeNum6  $xdim $ydim $zdim
             puts $meshFile "$nodeNum6  $xdim $ydim $zdim"
             set   nodeNum6   [expr $nodeNum6+1]
             
             node  $nodeNum6  $xdim $ydim_2 $zdim 
             puts $meshFile "$nodeNum6  $xdim $ydim_2 $zdim"
             set   nodeNum6   [expr $nodeNum6+1]
             

      }
   }
   
##################################################################################
###### Create dummy nodes to connect soil nodes to the abutment shell sodes ######
##################################################################################

# Left side
set   numnodebelow 1
set   incriment 0.1
set   L0    [expr $Lzone1+$Lzone2+$Lzone3+$WInfemb+$BP_emb/2.0]
set   dummynum  1000000
set   nodeNum301   $dummynum 

for {set i 1} {$i <= 4} {incr i 1} {
  for {set j 1} {$j <=  [expr $numYele_deck+1+4] } {incr j 1} {
      set m 1
      for {set k 1} {$k <= [expr $numZele4-$numnodebelow] } {incr k 1} {
	  
           if  { $i <= 2 } {
               if { $i == 1 } {
	           set xdim [expr $L0+$BP_emb/2.0]
		       } else {
		       if { $i == 2 } {
			   set xdim [expr $L0+$BP_emb/2.0+$WInfemb]
               }}		
		       }
			   
		   if  { $i == 3 } {
		       if  { $k ==  [expr $numZele4-$numnodebelow] } {
			       set xdim [expr $L0+$BP_emb/2.0+$WInfemb+$incriment]
				   } else {
		           set xdim [expr $L0+$BP_emb/2.0+$WInfemb+($Lzone5-$k*($zSize4*$sloperatio2))/2.0]
				   }
		        }
			   
		   if  { $i == 4 } {
		   		if  { $k ==  [expr $numZele4-$numnodebelow] } {
			        set xdim [expr $L0+$BP_emb/2.0+$WInfemb+2.0*$incriment]
				    } else {
		            set xdim [expr $L0+$BP_emb/2.0+$WInfemb+$Lzone5-$k*($zSize4*$sloperatio2)]
					}
		       }
		 
           if { $j <= 2	} { 	
		       if { $j == 1 } {
	           set ydim [expr $Bzone1+($zSize4/$sloperatio)*$k]
		       } else {
		       if { $j == 2 } {
			   set ydim [expr $Bzone1+($zSize4/$sloperatio)*$k+($Bzone2-($zSize4/$sloperatio)*$k)/2.0]
               }}	
			   }
		   if { $j > 2	&& $j <=  [expr $numYele_deck+1+2] } { 
		       set t [expr $j-2]
             if { $t <= [expr $numYele3+1] } {
             set ydim [expr $Bzone1+$Bzone2+(($t-1)*$ySize3)]
             } else {
             if { $t == [expr $numYele3+1+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
             } else {			 
	         if { $t > [expr $numYele3+1+1] && $t <= [expr $numYele3+1+1+2] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+($t-($numYele3+1+1))*$ySize4]
             } else {
	         if { $t == [expr ($numYele3+1+1+3)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
	         if { $t == [expr ($numYele3+1+1+4)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
            } else {
	      if { $t == [expr ($numYele3+1+1+5)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]	  
             } else {
	      if { $t == [expr ($numYele3+1+1+6)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]  	
             } else {			 
			 	############################	 
	  if { $t == [expr ($numYele3+9)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]  	
         } else {
	  if { $t == [expr ($numYele3+10)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]  	
         } else {		 
	  if { $t == [expr ($numYele3+11)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]  	
         } else {	
	  if { $t == [expr ($numYele3+12)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]  	
         } else {
	  if { $t == [expr ($numYele3+13)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]  		 
         } else {
	  if { $t == [expr ($numYele3+14)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {			 
	  if { $t == [expr ($numYele3+15)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {			 
	  if { $t == [expr ($numYele3+16)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]  	
         } else {	
	  if { $t == [expr ($numYele3+17)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {	
	  if { $t == [expr ($numYele3+18)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {	
	###################################	 
	
	      if { $t == [expr ($numYele3+7+$numXele11+1+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]	
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+2+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]		 
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+3+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]		
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+4+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	
             } else {
          if { $t == [expr ($numYele3+7+$numXele11+5+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+6+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
             } else {
	      if { $t == [expr $numYele3+7+$numXele11+7+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
             } else {
             set ydim [expr $B1+($t-($numYele3+7+$numXele11+7+1))*$ySize5]
             }}}}}}}}}}}}}}}}}}}}}}}}            
			 
             }			 

            if { $j >  [expr $numYele_deck+1+2] } { 
			    if { $j ==  [expr $numYele_deck+1+3] } {          
                   set ydim [expr $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k-($B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k-($B1+$Bzone5))/2.0]
			       } else {
			       set ydim   [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k]
				   }
			   }
			   
             set zdim   [expr $Czone1+$Czone2+$Czone3+$k*$zSize4]
			             
             node  $nodeNum301  $xdim $ydim $zdim
             puts $meshFile "$nodeNum301  $xdim $ydim $zdim"
             set   nodeNum301   [expr $nodeNum301+1]
             #set m [expr $m+1]	                             
          }	
	  
      }
   }    

### Right side
set numnodebelow 1
set incriment 0.1
set dummynum2  2000000
set   nodeNum302   $dummynum2 

for {set i 1} {$i <= 4} {incr i 1} {
  for {set j 1} {$j <=  [expr $numYele_deck+1+4] } {incr j 1} {
      set m 1
      for {set k 1} {$k <= [expr $numZele4-$numnodebelow] } {incr k 1} {
	  
           if  { $i <= 2 } {
               if { $i == 1 } {
	           set xdim [expr $L3-$BP_emb/2.0]
		       } else {
		       if { $i == 2 } {
			   set xdim [expr $L3-$BP_emb/2.0-$WInfemb]
               }}		
		       }
			   
		   if  { $i == 3 } {
		       if  { $k ==  [expr $numZele4-$numnodebelow] } {
			       set xdim [expr $L3-$BP_emb/2.0-$WInfemb-$incriment]
				   } else {
		           set xdim [expr $L3-$BP_emb/2.0-$WInfemb-($Lzone17-$k*($zSize4*$sloperatio2))/2.0]
				   }
		        }
			   
		   if  { $i == 4 } {
		   		if  { $k ==  [expr $numZele4-$numnodebelow] } {
			        set xdim [expr $L3-$BP_emb/2.0-$WInfemb-2.0*$incriment]
				    } else {
		            set xdim [expr $L3-$BP_emb/2.0-$WInfemb-$Lzone17+$k*($zSize4*$sloperatio2)]
					}
		       }
		
            
           if { $j <= 2	} { 	
		       if { $j == 1 } {
	           set ydim [expr $Bzone1+($zSize4/$sloperatio)*$k]
		       } else {
		       if { $j == 2 } {
			   set ydim [expr $Bzone1+($zSize4/$sloperatio)*$k+($Bzone2-($zSize4/$sloperatio)*$k)/2.0]
               }}	
			   }
		   if { $j > 2	&& $j <=  [expr $numYele_deck+1+2] } { 
		       set t [expr $j-2]


             if { $t <= [expr $numYele3+1] } {
             set ydim [expr $Bzone1+$Bzone2+(($t-1)*$ySize3)]
             } else {
             if { $t == [expr $numYele3+1+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb] 
             } else {			 
	         if { $t > [expr $numYele3+1+1] && $t <= [expr $numYele3+1+1+2] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+($t-($numYele3+1+1))*$ySize4]
             } else {
	         if { $t == [expr ($numYele3+1+1+3)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile]
             } else {
	         if { $t == [expr ($numYele3+1+1+4)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile/2.0]
            } else {
	      if { $t == [expr ($numYele3+1+1+5)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile]	  
             } else {
	      if { $t == [expr ($numYele3+1+1+6)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile]  	
             } else {
			 
			 	############################	 
	  if { $t == [expr ($numYele3+9)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11]  	
         } else {
	  if { $t == [expr ($numYele3+10)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$xSize10]  	
         } else {		 
	  if { $t == [expr ($numYele3+11)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile]  	
         } else {	
	  if { $t == [expr ($numYele3+12)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11]  	
         } else {
	  if { $t == [expr ($numYele3+13)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$xSize10]  		 
         } else {
	  if { $t == [expr ($numYele3+14)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {			 
	  if { $t == [expr ($numYele3+15)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {			 
	  if { $t == [expr ($numYele3+16)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$xSize10]  	
         } else {	
	  if { $t == [expr ($numYele3+17)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile]  	
         } else {	
	  if { $t == [expr ($numYele3+18)] } {
         set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11+$BP_pile+$xSize11]  	
         } else {	
	###################################	 
	
	      if { $t == [expr ($numYele3+7+$numXele11+1+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile]	
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+2+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile/2.0]		 
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+3+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile]		
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+4+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile]	
             } else {
          if { $t == [expr ($numYele3+7+$numXele11+5+1)] } {
              set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb/2.0]
             } else {
	      if { $t == [expr ($numYele3+7+$numXele11+6+1)] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb]
             } else {
	      if { $t == [expr $numYele3+7+$numXele11+7+1] } {
             set ydim [expr $Bzone1+$Bzone2+$Bzone3+$WInfemb+$BP_emb+$WInfpile+$BP_pile+$WInfpile+$Lzone11+$WInfpile+$BP_pile+$WInfpile+$BP_emb+$WInfemb]
             } else {
             set ydim [expr $B1+($t-($numYele3+7+$numXele11+7+1))*$ySize5]
             }}}}}}}}}}}}}}}}}}}}}}}}            
			 
             }			 

            if { $j >  [expr $numYele_deck+1+2] } { 
			    if { $j ==  [expr $numYele_deck+1+3] } {          
                   set ydim [expr $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k-($B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k-($B1+$Bzone5))/2.0]
			       } else {
			       set ydim   [expr  $B1+$Bzone5+$Bzone6-($zSize4/$sloperatio)*$k]
				   }
			   }
			   
             set zdim   [expr $Czone1+$Czone2+$Czone3+$k*$zSize4]
			             
             node  $nodeNum302  $xdim $ydim $zdim
             puts $meshFile "$nodeNum302  $xdim $ydim $zdim"
             set   nodeNum302   [expr $nodeNum302+1]
             #set m [expr $m+1]	                             
          }	
	  
      }
   }    

puts $meshFile "end coordinates"

	 
##################################
###### Create soil elements ######
##################################

puts $meshFile "Elements"

set m 0
set eleNum 0
for {set i 1} {$i <= $numXele} {incr i 1} {
  for {set j 1} {$j <= $numYele} {incr j 1} {
      for {set k 1} {$k <= $numZele} {incr k 1} {
	  
          set eleNum [expr $eleNum+1]
			  
          set n1  [expr ($i-1)*($numYnode*$numZnode)+$k+($j-1)*$numZnode] 
          set n2  [expr $n1+($numYnode*$numZnode)] 
          set n3  [expr $n2+$numZnode] 
          set n4  [expr $n1+$numZnode] 
          set n5  [expr $n1+1]
          set n6  [expr $n2+1]
          set n7  [expr $n3+1]
          set n8  [expr $n4+1]  	      
		
	 # User-Defined --- $k > ? values may change depending on soil layering

	 # Pile cap elements
     if { ($i >= $num_pile_X2 && $i <= [expr $num_pile_X3+3]) && ($j >= $num_pile_Y2 && $j <= $num_pile_Y4) && $k == 12 } {
          ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6 0.0 0.0 [expr $gravZ*$rhocap]
           element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6 0.0 0.0 [expr $gravZ*$rhocap]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6"
        } else {
     if { ($i >= [expr $num_pile_X2+3] && $i <= [expr $num_pile_X3]) && ($j > [expr $num_pile_Y2+2] && $j <= [expr $num_pile_Y4-3]) && $k == 13 } {
           ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6 0.0 0.0 [expr $gravZ*$rhocap]
           element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6 0.0 0.0 [expr $gravZ*$rhocap]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 6"
        } else {
		
	 # In situ soil material	
     if { $k > 11 } {
           ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1 0.0 0.0 [expr $gravZ*$rhos1]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1"
		   element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1 0.0 0.0 [expr $gravZ*$rhos1]
        } else {
     if { $k <= 11 &&  $k >= 7 } {
          ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 2 0.0 0.0 [expr $gravZ*$rhos2]
	       element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 2 0.0 0.0 [expr $gravZ*$rhos2]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 2"
        } else {
     if { $k < 7 &&  $k >=5 } {
           ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 3 0.0 0.0 [expr $gravZ*$rhos3]
	       element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 3 0.0 0.0 [expr $gravZ*$rhos3]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 3"
         } else {
     if { $k < 5 &&  $k > 2 } {
           ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 4 0.0 0.0 [expr $gravZ*$rhos4]
	       element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 4 0.0 0.0 [expr $gravZ*$rhos4]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 4"
         } else {
           ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 5 0.0 0.0 [expr $gravZ*$rhos5]
	       element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 5 0.0 0.0 [expr $gravZ*$rhos5]
           puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 5"
     }}}}}}
    }
  }
}

# embankment elements
for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele2 ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4 } {incr k 1} {
  
          set eleNum [expr $eleNum+1]
          if { $j == 1 } {     
             if { $k == 1 } {

                set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numYele1+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum+($i-1)*(2.0*$numZele4)+$k)]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($n6+1)]
                set n8  [expr int($n5+1)] 
                 	      
	            } else {

                set n1  [expr int($nodeNum+($i-1)*(2.0*$numZele4)+2*($k-1)-1)]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($n2+1)] 
                set n4  [expr int($n1+1)] 
                set n5  [expr int($n1+2)]
                set n6  [expr int($n2+2)]
                set n7  [expr int($n3+2)]
                set n8  [expr int($n4+2)]
                }
            } else {
               if { $k == 1 } {
                set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numYele1+2)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum+($i-1)*(2.0*$numZele4)+1+$k)]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($nodeNum1 + $i*(($numYele_deck+1)*$numZele4))]
                set n8  [expr int($nodeNum1+($i-1)*(($numYele_deck+1)*$numZele4))] 
	      
	            } else {
                
                set n1  [expr int($nodeNum+($i-1)*(2.0*$numZele4)+1+2*($k-1)-1)]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($nodeNum1 + $i*(($numYele_deck+1)*$numZele4)+($k-2))] 
                set n4  [expr int($nodeNum1 + ($i-1)*(($numYele_deck+1)*$numZele4)+($k-2))]  
                set n5  [expr int($n1+2)]
                set n6  [expr int($n2+2)]
                set n7  [expr int($n3+1)]
                set n8  [expr int($n4+1)]

                }
            }
	     ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
         element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

    }
  }
}
# to record outputs
set  stressstrainEMBL [expr $eleNum]

# embankment elements
for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele_deck ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4} {incr k 1} {
	  
          set eleNum [expr $eleNum+1]
          if { $k == 1 } {     
       
                set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numYele1+$numYele2+1)*$numZnode+($j-1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum1+($i-1)*($numYele_deck+1)*$numZele4+($j-1)*$numZele4)]
                set n6  [expr int($n5+(($numYele_deck+1)*$numZele4))]
                set n7  [expr int($n6+$numZele4)]
                set n8  [expr int($n5+$numZele4)] 
                 	      
          } else {
                
                set n1  [expr int($nodeNum1+($i-1)*($numYele_deck+1)*$numZele4+($k-2)+($j-1)*$numZele4)]
                set n2  [expr int($n1+(($numYele_deck+1)*$numZele4))]
                set n3  [expr int($n2+$numZele4)] 
                set n4  [expr int($n1+$numZele4)] 
                set n5  [expr $n1+1]
                set n6  [expr $n2+1]
                set n7  [expr $n3+1]
                set n8  [expr $n4+1]
          }

	     ##element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
         element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

    }
  }
}

# embankment elements
for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele6 ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4} {incr k 1} {
	  
          set eleNum [expr $eleNum+1]
          if { $j != 1 } {     
             if { $k == 1 } {
                set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numYele_dum+1+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum2+(2*$numZele4-1)+($i-1)*(2.0*$numZele4))]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($n6-1)]
                set n8  [expr int($n5-1)] 

                 	      
	            } else {
                
                set n1  [expr int(($nodeNum2-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4)-2*($k-2))]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($n2-1)] 
                set n4  [expr int($n1-1)] 
                set n5  [expr int($n1-2)]
                set n6  [expr int($n2-2)]
                set n7  [expr int($n3-2)]
                set n8  [expr int($n4-2)]

                }
            } else {
               if { $k == 1 } {
                set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numYele_dum+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int(($nodeNum1-1)+($i-1)*($numYele_deck+1)*$numZele4 + $numYele_deck*$numZele4+1)]                
                set n6  [expr int($n5+($numYele_deck+1)*$numZele4)]
                set n7  [expr int(($nodeNum2-1)+(2.0*$numZele4)+$i*(2.0*$numZele4))]
                set n8  [expr int(($nodeNum2-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4))]
            	      
	            } else {                
                set n1  [expr int(($nodeNum1-1) + $numYele_deck*$numZele4+($i-1)*(($numYele_deck+1)*$numZele4)+($k-1))]
                set n2  [expr int($n1+(($numYele_deck+1)*$numZele4))]
                set n3  [expr int(($nodeNum2-1)+(2.0*$numZele4)+$i*(2.0*$numZele4)-2*($k-2))]
                set n4  [expr int(($nodeNum2-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4)-2*($k-2))]  
                set n5  [expr int($n1+1)]
                set n6  [expr int($n2+1)]
                set n7  [expr int($n3-2)]
                set n8  [expr int($n4-2)]

                  }
            }

	    ## element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
		element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
        puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"
    }
  }
}


# embankment elements
set NodeNumTot [expr $numXELE_d*$numYnode*$numZnode]
   
 for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele2 ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4 } {incr k 1} {
  
          set eleNum [expr $eleNum+1]
          if { $j == 1 } {     
             if { $k == 1 } {

                set n1  [expr int($NodeNumTot+($i-1)*($numYnode*$numZnode)+($numYele1+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int(($nodeNum3-1)+($i-1)*(2.0*$numZele4)+$k)]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($n6+1)]
                set n8  [expr int($n5+1)] 
                 	      
	            } else {

                set n1  [expr int(($nodeNum3-1)+($i-1)*(2.0*$numZele4)+2*($k-1)-1)]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($n2+1)] 
                set n4  [expr int($n1+1)] 
                set n5  [expr int($n1+2)]
                set n6  [expr int($n2+2)]
                set n7  [expr int($n3+2)]
                set n8  [expr int($n4+2)]
                }
            } else {
               if { $k == 1 } {
                set n1  [expr int($NodeNumTot+($i-1)*($numYnode*$numZnode)+($numYele1+2)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int(($nodeNum3-1)+($i-1)*(2.0*$numZele4)+1+$k)]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($nodeNum4 + $i*(($numYele_deck+1)*$numZele4))]
                set n8  [expr int($nodeNum4+($i-1)*(($numYele_deck+1)*$numZele4))] 
	      
	            } else {
                
                set n1  [expr int(($nodeNum3-1)+($i-1)*(2.0*$numZele4)+1+2*($k-1)-1)]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($nodeNum4 + $i*(($numYele_deck+1)*$numZele4)+($k-2))] 
                set n4  [expr int($nodeNum4 + ($i-1)*(($numYele_deck+1)*$numZele4)+($k-2))]  
                set n5  [expr int($n1+2)]
                set n6  [expr int($n2+2)]
                set n7  [expr int($n3+1)]
                set n8  [expr int($n4+1)]

                }
            }
	     ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
         element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

    }
  }
}

set  stressstrainEMBR [expr $eleNum]

# embankment elements
for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele_deck ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4} {incr k 1} {
	  
          set eleNum [expr $eleNum+1]
          if { $k == 1 } {     
       
                set n1  [expr int($NodeNumTot+($i-1)*($numYnode*$numZnode)+($numYele1+$numYele2+1)*$numZnode+($j-1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum4+($i-1)*($numYele_deck+1)*$numZele4+($j-1)*$numZele4)]
                set n6  [expr int($n5+(($numYele_deck+1)*$numZele4))]
                set n7  [expr int($n6+$numZele4)]
                set n8  [expr int($n5+$numZele4)] 
                	      
          } else {
                
                set n1  [expr int($nodeNum4+($i-1)*($numYele_deck+1)*$numZele4+($k-2)+($j-1)*$numZele4)]
                set n2  [expr int($n1+(($numYele_deck+1)*$numZele4))]
                set n3  [expr int($n2+$numZele4)] 
                set n4  [expr int($n1+$numZele4)] 
                set n5  [expr $n1+1]
                set n6  [expr $n2+1]
                set n7  [expr $n3+1]
                set n8  [expr $n4+1]
          }

	    ### element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
		element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
        puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

    }
  }
}


# embankment elements
for {set i 1} {$i <= $numXele_emb} {incr i 1} {
  for {set j 1} {$j <= [expr $numYele6 ]} {incr j 1} {
      for {set k 1} {$k <= $numZele4} {incr k 1} {
	  
          set eleNum [expr $eleNum+1]
          if { $j != 1 } {     
             if { $k == 1 } {
                set n1  [expr int($NodeNumTot+($i-1)*($numYnode*$numZnode)+($numYele_dum+1+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int($nodeNum5+(2*$numZele4-1)+($i-1)*(2.0*$numZele4))]
                set n6  [expr int($n5+(2.0*$numZele4))]
                set n7  [expr int($n6-1)]
                set n8  [expr int($n5-1)] 
                 	      
	            } else {
                
                set n1  [expr int(($nodeNum5-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4)-2*($k-2))]
                set n2  [expr int($n1+(2.0*$numZele4))]
                set n3  [expr int($n2-1)] 
                set n4  [expr int($n1-1)] 
                set n5  [expr int($n1-2)]
                set n6  [expr int($n2-2)]
                set n7  [expr int($n3-2)]
                set n8  [expr int($n4-2)]
                }
            } else {
               if { $k == 1 } {
                set n1  [expr int($NodeNumTot+($i-1)*($numYnode*$numZnode)+($numYele_dum+1)*$numZnode)] 
                set n2  [expr int($n1+($numYnode*$numZnode))] 
                set n3  [expr int($n2+$numZnode)] 
                set n4  [expr int($n1+$numZnode)] 
                set n5  [expr int(($nodeNum4-1)+($i-1)*($numYele_deck+1)*$numZele4 + $numYele_deck*$numZele4+1)]                
                set n6  [expr int($n5+($numYele_deck+1)*$numZele4)]
                set n7  [expr int(($nodeNum5-1)+(2.0*$numZele4)+$i*(2.0*$numZele4))]
                set n8  [expr int(($nodeNum5-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4))]
                 	      
	            } else {                
                set n1  [expr int(($nodeNum4-1) + $numYele_deck*$numZele4+($i-1)*(($numYele_deck+1)*$numZele4)+($k-1))]
                set n2  [expr int($n1+(($numYele_deck+1)*$numZele4))]
                set n3  [expr int(($nodeNum5-1)+(2.0*$numZele4)+$i*(2.0*$numZele4)-2*($k-2))]
                set n4  [expr int(($nodeNum5-1)+(2.0*$numZele4)+($i-1)*(2.0*$numZele4)-2*($k-2))]  
                set n5  [expr int($n1+1)]
                set n6  [expr int($n2+1)]
                set n7  [expr int($n3-2)]
                set n8  [expr int($n4-2)]

                  }
            }
	     ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
         element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]			 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

    }
  }
}

# embankment elements
 for {set i 1} { $i <= 4 } {incr i 1} {
   for {set j 1} {$j <= [expr $numYele_deck+4] } {incr j 1} {
       for {set k 1} {$k <=  [expr $numZele4-$numnodebelow]   } {incr k 1} {

           set eleNum [expr $eleNum+1]
           if { $i == 1 && $j == 1 && $k ==1 } {     

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int(($nodeNum1-1)-2.0*$numZele4+1+2.0*($k-1))]
                 set n6  [expr int($dummynum)]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int($n5+1)] 
              }
           if { $i == 1 && $j == 2 && $k ==1 } {     

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int(($nodeNum1-1)-2.0*$numZele4+1+2.0*($k-1)+1)]
                 set n6  [expr int($dummynum+($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1))]  
              }  

           if { $i == 1 && $k ==1 } {   
               if { $j > 2 && $j <= [expr $numYele_deck+2] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1)+($j-3)*$numZele4)]
                 set n6  [expr int($dummynum+($j-1)*($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int($n5+$numZele4)]  
              }        
			   if { $j == [expr $numYele_deck+3] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1)+($j-3)*$numZele4)]
                 set n6  [expr int($dummynum+($j-1)*($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($nodeNum3-1))]   
              }   
			   if { $j == [expr $numYele_deck+4] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int(($nodeNum3-1))]
                 set n6  [expr int($dummynum+($j-1)*($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($n5-1))]   
              }   
             }		

           if { $i == 1 && $j == 1 && $k >1 } {     

                 set n1  [expr int(($nodeNum1-1)-2.0*$numZele4+1+2.0*($k-2))]
                 set n2  [expr int($dummynum+($k-2))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))]
                 set n4  [expr int($n5+1)] 
                 set n5  [expr int($n1+2)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int($n5+1)] 
              }
           if { $i == 1 && $j == 2 && $k >1 } {     

                 set n1  [expr int(($nodeNum1-1)-2.0*$numZele4+1+2.0*($k-2)+1)]
                 set n2  [expr int($dummynum+($k-2)+($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))]
                 set n4  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1 +($k-2)))] 
                 set n5  [expr int($n1+2)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($n4+1))]  
              }  

           if { $i == 1 && $k >1 } {   
               if { $j > 2 && $j <= [expr $numYele_deck+2] } { 		   

                 set n1  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1)+($k-2)+($j-3)*$numZele4)]
                 set n2  [expr int($dummynum+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))] 
                 set n4  [expr int($n1+$numZele4)] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n3+1)]
                 set n8  [expr int($n4+1)]  
              }        
			   if { $j == [expr $numYele_deck+3] } { 		   

                 set n1  [expr int(($nodeNum2-1-($numYele_deck+1)*$numZele4+1)+($k-2)+($j-3)*$numZele4)]
                 set n2  [expr int($dummynum+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))] 
                 set n4  [expr int($nodeNum3-1-2*($k-2))] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($n4-2))]   
              }   
			   if { $j == [expr $numYele_deck+4] } { 		   

                 set n1  [expr int($nodeNum3-1-2*($k-2))]
                 set n2  [expr int($dummynum+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))]  
                 set n4  [expr int($nodeNum3-1-2*($k-2)-1)] 
                 set n5  [expr int(($n1-2))]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int(($n5-1))]   
              }   
             }			  
			 
			 if { $i > 1 && $k==1 } {    

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXele1+$numXele2+$numXele3+2)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int($dummynum+($i-2)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($j-1)*($numZele4-$numnodebelow))]
                 set n6  [expr int($n5+($numYele_deck+4+1)*($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]  
               }  			 
			  if { $i > 1 && $k > 1 } {    

                 set n1  [expr int($dummynum+($i-2)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n2  [expr int($n1+($numYele_deck+4+1)*($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n3+1)]
                 set n8  [expr int($n4+1)] 
               }  			 
			 
	     ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
         element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

     }
   }
 }

 
# embankment elements
 for {set i 1} { $i <= 4 } {incr i 1} {
   for {set j 1} {$j <= [expr $numYele_deck+4] } {incr j 1} {
       for {set k 1} {$k <=  [expr $numZele4-$numnodebelow]   } {incr k 1} {

           set eleNum [expr $eleNum+1]
           if { $i == 1 && $j == 1 && $k ==1 } {     

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
				 set n5  [expr int($dummynum2)]
                 set n6  [expr int(($nodeNum3-1)+1+2.0*($k-1))]
				 set n7  [expr int($n6+1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]
 
              }
           if { $i == 1 && $j == 2 && $k ==1 } {     

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
				 set n5  [expr int($dummynum2+($numZele4-$numnodebelow))]
                 set n6  [expr int(($nodeNum3-1)+1+2.0*($k-1)+1)]
				 set n7  [expr int(($nodeNum4-1)+1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]
              }  

           if { $i == 1 && $k ==1 } {   
               if { $j > 2 && $j <= [expr $numYele_deck+2] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
				 set n5  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow))]
				 set n6  [expr int(($nodeNum4-1)+1+($j-3)*$numZele4)]
                 set n7  [expr int($n6+$numZele4)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]  
                 }        
			   if { $j == [expr $numYele_deck+3] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow))]
                 set n6  [expr int(($nodeNum4-1)+1+$numZele4*($numYele_deck))]
                 set n7  [expr int($nodeNum5-1+2*$numZele4)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]    
              }   
			   if { $j == [expr $numYele_deck+4] } { 		   

                 set n1  [expr int(($i-1)*($numYnode*$numZnode)+($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow))]
                 set n6  [expr int($nodeNum5-1+2*$numZele4)]
                 set n7  [expr int($n6-1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]     
              }   
              }		
		 
           if { $i == 1 && $j == 1 && $k >1 } {     

                 set n1  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow)+($k-2))]
                 set n2  [expr int(($nodeNum3-1)+1+2.0*($k-2))]
                 set n3  [expr int($n2+2*($numZele4))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+2)]
                 set n7  [expr int($n6+1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))] 
              }
           if { $i == 1 && $j == 2 && $k >1 } {     

                 set n1  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow)+($k-2))]
                 set n2  [expr int(($nodeNum3-1)+1+2.0*($k-2)+1)]
                 set n3  [expr int(($nodeNum4-1)+1+($k-2))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))]  
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+2)]
                 set n7  [expr int($n3+1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]   
              }  

           if { $i == 1 && $k >1 } {   
               if { $j > 2 && $j <= [expr $numYele_deck+2] } { 		   

                 set n1  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow)+($k-2))]
                 set n2  [expr int(($nodeNum4-1)+1+($k-2)+($j-3)*$numZele4)]
                 set n3  [expr int($n2+($numZele4))] 
                 set n4  [expr int($n1+($numZele4-$numnodebelow))] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n3+1)]
                 set n8  [expr int($n4+1)]  
              }        
			   if { $j == [expr $numYele_deck+3] } { 		   

                 set n1  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow)+($k-2))]
                 set n2  [expr int(($nodeNum4-1)+1+($k-2)+($j-3)*$numZele4)]
                 set n3  [expr int($nodeNum5-1+2.0*$numZele4-2.0*($k-2))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))]  
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n3-2)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))] 
              }   
			   if { $j == [expr $numYele_deck+4] } { 		   

                 set n1  [expr int($dummynum2+($j-1)*($numZele4-$numnodebelow)+($k-2))]
                 set n2  [expr int(($nodeNum5-1)+2.0*$numZele4-2.0*($k-2))]
                 set n3  [expr int(($n2-1))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))]  
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2-2)]
                 set n7  [expr int($n6-1)]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]  
              }  
             }			  
			 
			 if { $i > 1 && $k==1 } {    

                 set n1  [expr int(($numXELE_d-$i)*($numYele+1)*$numZnode+($numYele1+$j)*$numZnode)] 
                 set n2  [expr int($n1+($numYnode*$numZnode))] 
                 set n3  [expr int($n2+$numZnode)] 
                 set n4  [expr int($n1+$numZnode)] 
                 set n5  [expr int($dummynum2+($i-1)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($j-1)*($numZele4-$numnodebelow))]
                 set n6  [expr int($dummynum2+($i-2)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($j-1)*($numZele4-$numnodebelow))]
                 set n7  [expr int($n6+($numZele4-$numnodebelow))]
                 set n8  [expr int($n5+($numZele4-$numnodebelow))]  
               }  			 
			  if { $i > 1 && $k > 1 } {    

                 set n1  [expr int($dummynum2+($i-1)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n2  [expr int($dummynum2+($i-2)*($numYele_deck+4+1)*($numZele4-$numnodebelow)+($k-2)+($j-1)*($numZele4-$numnodebelow))]
                 set n3  [expr int($n2+($numZele4-$numnodebelow))]
                 set n4  [expr int($n1+($numZele4-$numnodebelow))] 
                 set n5  [expr int($n1+1)]
                 set n6  [expr int($n2+1)]
                 set n7  [expr int($n3+1)]
                 set n8  [expr int($n4+1)] 
               }  			 
			 
	     ###element stdBrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		
		 element SSPbrick $eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7 0.0 0.0 [expr $gravZ*$rhosemb]		 
         puts $meshFile "$eleNum $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 7"

     }
   }
 }
 
 
puts $meshFile "end elements"


##################################
###### Create soil elements ######
##################################

# base of the model
for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {

       fix [expr 1+($i-1)*$numYnode*$numZnode+($j-1)*$numZnode] 1 1 1 
         
      }
    }


# Lateral Boundaries 
for {set k 2 } { $k <= $numZnode } { incr k 1 } {
  for {set j 2} {$j <=  $numYnode } { incr j 1 } {

      equalDOF $k [expr ($numXnode-1)*($numZnode*$numYnode)+$k+(($j-2)*$numZnode)]  1  2  3      ;# right side yz boundary 
      equalDOF $k [expr ($j-1)*$numZnode+$k]  1   2  3                                           ;# left  side yz boundary                                                                                
        }

      equalDOF $k  [expr int($numXnode*$numZnode*$numYnode-$numZnode+$k)]  1 2 3

  for {set i 2} {$i <=  [expr ($numXnode-1)]} {incr i 1} {
     equalDOF $k [expr ($i-1)*($numZnode*$numYnode)+$k]      1    2  3                            ;# right side xz boundary
     equalDOF $k [expr $numZnode*($numYnode-1)+($i-1)*($numZnode*$numYnode)+$k]      1    2  3    ;# right side xz boundary
       }
}


# Lateral Boundaries for embankments
# BCs at X=0 and X= model length
set  nodeNum_emb_BC  [expr $nodeNum+1]
set m 0

for {set k 1 } { $k <= [expr $numZele4] } { incr k 1 } {

      equalDOF    $nodeNum_emb_BC     [expr $nodeNum_emb_BC+1]        1  2  3 
      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum2-1)+2.0*$numZele4-1-$m)]                    1 2 3   
      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum2-1)+2.0*$numZele4-1-$m+1)]                  1 2 3   
  
      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum3-1)+2.0*$numZele4*($numXnode_emb-1)+1+$m)]                1 2 3
      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum3-1)+2.0*$numZele4*($numXnode_emb-1)+1+$m+1)]              1 2 3

      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum5-1)+2.0*$numZele4*($numXnode_emb)-1-$m)]                  1 2 3
      equalDOF    $nodeNum_emb_BC     [expr int(($nodeNum5-1)+2.0*$numZele4*($numXnode_emb)-1-$m+1)]                1 2 3

      for {set j 1} {$j <=  [expr $numYele_deck+1] } { incr j 1 } {

            equalDOF    $nodeNum_emb_BC     [expr int($nodeNum1+($j-1)*$numZele4+($k-1))]   1 2 3                        
            equalDOF    $nodeNum_emb_BC     [expr int($nodeNum4+($numXnode_emb-1)*($numYele_deck+1)*$numZele4+($j-1)*$numZele4+($k-1))]   1 2  3
        
            }

      set    nodeNum_emb_BC    [expr $nodeNum_emb_BC+2]   
      set    m   [expr $m+2]  
     }