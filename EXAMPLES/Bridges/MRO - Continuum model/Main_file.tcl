
# Created by Amin Rahmani - April 3, 2013
# Copyright © 2013 Amin Rahmani. All rights reserved.

#========================================================================================
#                               ****   Notes  ****
# This code performs a 3D seismic analysis for a full-scale two-span bridge, i.e. Meloland Road Overpass, in El Centro, CA.
# You need to make changes on the lines that the comment "User-Defined" appears. Please go through all the command lines before running it.
# Units are in Metric (KN, ton, sec, m). Check the units if you work with Imperial units. 
# When rigid boundary condition is used, the input motions should be available in two separate .txt files. Time in "SCTime.txt" and Displacement in "Disp_X,Y.txt".
# When compliant boundary condition is used, the input motions should be available in two separate .txt files. Time in "Time.txt" and Velocity in "Velocity_X,Y.txt".
# Outputs will be saved in a folder titled "Outputs".
#========================================================================================

wipe all

#after 180000

# Measures the total analysis execution time
set startT [clock seconds]

# Creates input file for GiD software
set meshFile [open soilmesh.flavia.msh w]

# Load steps for analysis
set NumIncr1  10
set NumIncr2  100

################################
###### Create Soil Domain ######
################################

model basic -ndm 3 -ndf 3
source Soil_Domain.tcl

########################################################
###### Geostatic Analysis - Soil in Elastic State ######
########################################################
 
 #integrator LoadControl [expr 1.0/$NumIncr1] 
 integrator LoadControl 1 1 1 1
 numberer Plain
 constraints Penalty 1.e12 1.e12
 test NormDispIncr 1.0e-4 50 1
 algorithm KrylovNewton
#system SparseGeneral
 # system Mumps
 system SparseSYM
 analysis Static

 for {set numIncr 1} {$numIncr <= $NumIncr1 } {incr numIncr 1} {
       puts "##### Elastic Geostatic Analysis step : $numIncr #####";
       analyze 1 
     }

##########################################
###### Switch soil to plastic state ######
##########################################

 updateMaterialStage -material 1 -stage 1  ; #Use-defined --- Specify the material numbers
 updateMaterialStage -material 2 -stage 1  ; #Use-defined --- Specify the material numbers
 updateMaterialStage -material 3 -stage 1  ; #Use-defined --- Specify the material numbers
 updateMaterialStage -material 4 -stage 1  ; #Use-defined --- Specify the material numbers
 updateMaterialStage -material 5 -stage 1  ; #Use-defined --- Specify the material numbers
 updateMaterialStage -material 7 -stage 1  ; #Use-defined --- Specify the material numbers

 for {set numIncr 1} {$numIncr <= $NumIncr1 } {incr numIncr 1} {
       puts "##### Plastic Geostatic Analysis step : $numIncr #####";
       analyze 1 
     }

############################
###### Reset Analyses ######
############################
 wipeAnalysis
 remove recorders
 loadConst -time 0.0

###############################################################################
###### Define Pile Group at the Pier - Sections and Material Properties  ######
###############################################################################

model basic -ndm 3 -ndf 6

# properties of piles
set D    0.32  					 ; #Use-defined --- Pile diameter
set t    0.00  					 ; #Use-defined --- Pile wall thickness
set A    [expr 3.1416*$D*$D/4.0] ; #Use-defined --- Pile cross-sectional area
set E    1.24e+7				 ; #Use-defined --- Pile Young's modulus 				 
set Fy   25.5e+4				 ; #Use-defined --- Pile yield strength
set G    [expr $E/2.0/(1+0.20)]  ; #Use-defined --- Pile Material Shear Modulus
set rhop 0.2					 ; #Use-defined --- Pile Material unit weight
set massDens_pilepier [expr $rhop*$A]

# properties of rigid connections horizontally connecting pile nodes to the surrounding soils
set Acon    [expr (3.1416*$D*$D/4.0)] 				  ; # Rigid elements
set Econ    1.24e+7 				  				  ; # Rigid elements  
set Fycon   25.5e+4 				 				  ; # Rigid elements 
set Gcon    [expr $Econ/2.0/(1+0.20)] 				  ; # Rigid elements 
set Izcon   [expr 1000.0*0.25*3.1416*(pow([expr $D/2.] , 4))] 				  ; # Rigid elements 
set Iycon   [expr 1000.0*0.25*3.1416*(pow([expr $D/2.] , 4))] 				  ; # Rigid elements 
set Jcon    [expr $Izcon+$Iycon]
uniaxialMaterial Elastic 7000 $E 

section Fiber 1 { 
  patch circ  7000 10 10 0.0 0.0 0.0 [expr $D/2] 0.0 360.0 ; #Use-defined --- Pile Cross Section
  }

source Piles_for_pier.tcl


###################################################################################
###### Define Pile Group at the abutment - Sections and Material Properties  ######
###################################################################################

set D_abut    0.32 								; #Use-defined --- Pile diameter
set rhop_abut 0.2  								; #Use-defined --- Pile Material unit weight
set A_abut    [expr 3.1416*$D_abut*$D_abut/4.0] ; #Use-defined --- Pile cross-sectional area
set massDens_pileabut [expr $rhop_abut*$A_abut]

source Piles_for_abutment.tcl

############################################################
###### Define Pier - Section and Material Properties  ######
############################################################

set Hpier     6.0 							 ; #Use-defined --- Pier height
set psize     1.0 							 ; #Use-defined --- element size
set pierD     1.52 							 ; #Use-defined --- Pier diameter
set Epier     20.0e+6						 ; #Use-defined --- concrete Young's modulus 
set Ebar      200.0e+6						 ; #Use-defined --- bar Young's modulus 
set Apier     [expr 0.25*3.1416*pow($pierD,2)] ; #Use-defined --- pier cross-sectional area 
set rhopier   2.45							   ; #Use-defined --- pier material density
set Fypier    455.0e+3						   ; #Use-defined --- pier yield strength
set areaFiber 0.002581						   ; #Use-defined 
set massDens_pier [expr $rhopier*$Apier]

set fc -34473.8							       ; #Use-defined --- concrete strength
set epsc0  [expr 2.0*$fc/$Epier]


uniaxialMaterial Elastic 8000 $Epier ; #Use-defined
uniaxialMaterial Elastic 8005 $Ebar; #Use-defined
uniaxialMaterial Steel01 8001 $Fypier $Ebar 0.0; #Use-defined
uniaxialMaterial Steel02 8002 $Fy $Ebar 0.0 18 0.925 0.15; #Use-defined 
uniaxialMaterial Concrete01  8003  -34473.8  -0.004   -21000.0  -0.014 ; #Use-defined
uniaxialMaterial Concrete01  8013  -27600.0  -0.002   0.0  -0.008; #Use-defined

section Fiber 2 { 
  patch circ  8003 50 50 0.0 0.0 0. [expr $pierD/2.0-0.074] 0.0 360.0; #Use-defined
  patch circ  8013 50 50 0.0 0.0 [expr $pierD/2.0-0.074] [expr $pierD/2.0] 0.0 360.0; #Use-defined
  layer circ  8001 18 $areaFiber 0.0 0.0 0.686 0.0 360.0 ; #Use-defined
  }

source Pier.tcl


################################################################
###### Define the Deck - Section and Material Properties  ######
################################################################

set  Ldeck      [expr 2.0*$BP_emb/2.0+2.0*$WInfemb+$Lzone5+$Lzone6+$Lzone7+$Lzone8+$Lzone9+2.0*$WInfpile+2.0*$BP_pile+2.0*$WInfpile+$Lzone11+$Lzone13+$Lzone14+$Lzone15+$Lzone16+$Lzone17] ; #Total length of the deck 
set  Wdeck      [expr $Bzone3+2.0*$WInfemb+2.0*$BP_emb+4.0*$WInfpile+2.0*$BP_pile+$Lzone11+$Bzone5] ; #Total width of the deck
set Hdeck       1.25 ; #Use-defined --- Deck Thickness
puts "### Length of Deck = $Ldeck ####"
puts "### Width of Deck = $Wdeck ####"
set Adeck      [expr $Wdeck*$Hdeck]  
set Edeck      20.0e+6				 ; #Use-defined --- Deck Young's modulus
set noudeck    0.2					 ; #Use-defined --- Deck material Poisson's ratio
set rhodeck    2.450				 ; #Use-defined --- Deck material density

set Dzone1 [expr $Ldeck/2.0-int($Ldeck/2.0)]
set Dzone3 $Dzone1
set Dzone2 [expr $Ldeck-$Dzone1-$Dzone3]

set xSize_D1  $Dzone1
set xSize_D3  $Dzone3
set xSize_D2  1.0 ; #Use-defined --- Deck element size

set   numXele_D    [expr int($Dzone1/$xSize_D1+$Dzone2/$xSize_D2+$Dzone3/$xSize_D3)]
set   numYele_D    [expr int($numYele3+7+$numXele11+7+$numYele5)]

############## deck x section
nDMaterial ElasticIsotropic 3999 $Edeck $noudeck $rhodeck ; #Use-defined
section PlateFiber 32000 3999 $Hdeck 					  ; #Use-defined

source Deck.tcl

#####################################################################
###### Define the abutments - Section and Material Properties  ######
#####################################################################
set  abutHeight   [expr $Czone4-$Df_abut] ; #abutment wall height
set  abutWidth    $Wdeck 				  ; #abutment wall width
set  abutthick    0.46 					  ; #abutment wall thickness
set  abutzsize  $zSize4
set  numXeleabut 1.0
set  numYeleabut $numYele_D
set  numZeleabut [expr $abutHeight/$abutzsize]
set  Eabut        24.80e+6				 ; #User-defined --- abutment wall Young's modulus
set  nouabut      0.2					 ; #User-defined --- abutment wall Poisson's ratio
set  rhoabut      2.45					 ; #User-defined --- abutment wall mass density
set abutthick_wing 0.305

nDMaterial ElasticIsotropic 4999 $Eabut $nouabut $rhoabut ; #Use-defined
section PlateFiber 44000 4999 $abutthick 				  ; #Use-defined
nDMaterial ElasticIsotropic 5999 $Eabut $nouabut $rhoabut ; #Use-defined
section PlateFiber 54000 5999 $abutthick_wing			  ; #Use-defined

source Abutments.tcl


##### Close the input file for GID software
close $meshFile
###########################################

#############################################################################
###### Apply deck tributary weight on the pier and the abutment walls  ######
#############################################################################

pattern Plain 500 "Linear" {
  load [ expr int($nodeNum22+$Hpier/$psize+$Df/$zSize3)] 0.0 0.0 -9614.0 0.0 0.0 0.0 ; #Use-defined
  }
  
pattern Plain 550 "Linear" {

   for { set i 1 } { $i <=  [expr $numYele_D+1]  } { incr i 1 } {
     set   nodedeck_L [expr int($DeckLnode-1+$i)]
     set   nodedeck_R [expr int($DeckRnode-1+$i)]
	 load $nodedeck_L  0.0 0.0 [expr -9614.0/2.0/($numYele_D+1)] 0.0 0.0 0.0 ; #Use-defined
	 load $nodedeck_R  0.0 0.0 [expr -9614.0/2.0/($numYele_D+1)] 0.0 0.0 0.0 ; #Use-defined

	   }
}

########################################################
###### Define output files  for Gravaity Analysis ######
########################################################

set dataDir bridgedata
file mkdir $dataDir

set  numElE_b    [expr $numXele1+$numXele2+$numXele3+4+$numXele5+$numXele6+$numXele7]

eval "recorder Element -file $dataDir/stress1GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strain1GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stress2GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strain2GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stress3GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strain3GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+$num_pile_Y1*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stress4GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strain4GRAV.out   -time     -eleRange   [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numElE_b+$numXele8+$numXele9+4+$numXele11+3)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stressFFGRAV.out   -time     -eleRange   [expr int(($numELE_c)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numELE_c)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strainFFGRAV.out   -time     -eleRange   [expr int(($numELE_c)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+1)]    [expr int(($numELE_c)*($numYnode-1)+($num_pile_Y1+12)*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stressFF_2GRAV.out   -time     -eleRange   [expr int(($numELE_c)*($numYnode-1)+(2.0)*($numZnode-1)+1)]    [expr int(($numELE_c)*($numYnode-1)+(2.0)*($numZnode-1)+($numZnode-1))]     stress"
eval "recorder Element -file $dataDir/strainFF_2GRAV.out   -time     -eleRange   [expr int(($numELE_c)*($numYnode-1)+(2.0)*($numZnode-1)+1)]    [expr int(($numELE_c)*($numYnode-1)+(2.0)*($numZnode-1)+($numZnode-1))]     strain"

eval "recorder Element -file $dataDir/stressEMB_LGRAV.out   -time     -eleRange     [expr  int($stressstrainEMBL+($numXele_emb-1)*$numYele_deck*$numZele4+($numYele_deck/2.0-1)*$numZele4+1)]  [expr  int($stressstrainEMBL+($numXele_emb-1)*$numYele_deck*$numZele4+($numYele_deck/2.0-1)*$numZele4+$numZele4)]   stress"
eval "recorder Element -file $dataDir/strainEMB_LGRAV.out   -time     -eleRange    [expr  int($stressstrainEMBL+($numXele_emb-1)*$numYele_deck*$numZele4+($numYele_deck/2.0-1)*$numZele4+1)]  [expr  int($stressstrainEMBL+($numXele_emb-1)*$numYele_deck*$numZele4+($numYele_deck/2.0-1)*$numZele4+$numZele4)]    strain"

eval "recorder Element -file $dataDir/stressEMB_RGRAV.out   -time     -eleRange  [expr  int($stressstrainEMBR+($numYele_deck/2.0-1)*$numZele4+1)]  [expr  int($stressstrainEMBR+($numYele_deck/2.0-1)*$numZele4+$numZele4)] stress"
eval "recorder Element -file $dataDir/strainEMB_RGRAV.out   -time     -eleRange  [expr  int($stressstrainEMBR+($numYele_deck/2.0-1)*$numZele4+1)]  [expr  int($stressstrainEMBR+($numYele_deck/2.0-1)*$numZele4+$numZele4)]      strain"

###############################################################
###### Apply deck weight and reach the equilibrium state ######
###############################################################

 integrator LoadControl [expr 1.0/$NumIncr2] 
 #integrator LoadControl 1 1 1 1
 numberer Plain
 constraints Penalty 1.e12 1.e12
 test NormDispIncr 1.0e-4 50 1
 algorithm KrylovNewton
 system SparseGeneral
 #system Mumps
 system SparseSYM
 analysis Static

 for {set numIncr 1} {$numIncr <= $NumIncr2 } {incr numIncr 1} {
       puts "##### System Equilibrium Analysis step : $numIncr #####";
       analyze 1 
     }
	 
######################################
###### Perform Seismic Analysis ######
######################################

wipeAnalysis
remove recorders
loadConst -time 0.0

############################################################################################
###### Remove Constraints at the model base in the directions that GM will be applied ######
############################################################################################

for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {

        remove sp [expr 1+($i-1)*$numYnode*$numZnode+($j-1)*$numZnode]  1
        remove sp [expr 1+($i-1)*$numYnode*$numZnode+($j-1)*$numZnode]  2 
         
      }
    }

#######################################################
###### Define output files  for Seismic Analysis ######
#######################################################

set  pilezone    [expr $numELE_b+$numXele8+$numXele9+4+$numXele11/2]
set  pilezone2   [expr $numELE_b+$numXele8+$numXele9+4+$numXele11+4]
set  numElE_b    [expr $numXele1+$numXele2+$numXele3+4+$numXele5+$numXele6+$numXele7]

# Records acceleration time history at the free-field of embankment - Left side
set nodeList_emb_L2 {}
for {set j 1} { $j <= [expr $numYeleabut+1]} {incr j 1} {
      lappend  nodeList_emb_L2  [expr int($nodeNum1-1+1.0*($numYele_deck+1)*$numZele4+$numZele4+($j-1)*$numZele4)]
    }
eval "recorder Node -file $dataDir/Emb_L_accel_FF.out        -time    -node  $nodeList_emb_L2  -dof 1 2 accel"


# Records acceleration time history at the free-field of embankment - Right side
set nodeList_emb_R2 {}
for {set j 1} { $j <= [expr $numYeleabut+1]} {incr j 1} {
     lappend  nodeList_emb_R2  [expr int($nodeNum4-1+$numXnode_emb*($numYele_deck+1)*$numZele4-2.0*($numYele_deck+1)*$numZele4+$numZele4+($j-1)*$numZele4)]
    }
eval "recorder Node -file $dataDir/Emb_R_accel_FF.out        -time    -node  $nodeList_emb_R2  -dof 1 2 accel"

# Records acceleration, displacement, and reaction force time histories at the interface of deck and abutmenet wall - Left abutment
set nodeList_deck_L {}
for {set j 1} { $j <= [expr $numYeleabut+1]} {incr j 1} {
     lappend  nodeList_deck_L  [expr int($DeckLnode+($j-1))]
    }

eval "recorder Node -file $dataDir/Deck_L_accel.out        -time    -node  $nodeList_deck_L  -dof 1 2 accel"	
eval "recorder Node -file $dataDir/Deck_L_Disp.out         -time    -node  $nodeList_deck_L  -dof 1 2 disp"
eval "recorder Node -file $dataDir/Deck_L_reaction.out     -time    -node  $nodeList_deck_L  -dof 1 2 reaction"

# Records acceleration, displacement, and reaction force time histories at the interface of deck and abutmenet wall - Right abutment
set nodeList_deck_R {}
for {set j 1} { $j <= [expr $numYeleabut+1]} {incr j 1} {
     lappend  nodeList_deck_R  [expr int($DeckRnode+($j-1))]
    }

eval "recorder Node -file $dataDir/Deck_R_accel.out        -time    -node  $nodeList_deck_R  -dof 1 2 accel"
eval "recorder Node -file $dataDir/Deck_R_Disp.out         -time    -node  $nodeList_deck_R  -dof 1 2 disp"
eval "recorder Node -file $dataDir/Deck_R_reaction.out     -time    -node  $nodeList_deck_R  -dof 1 2 reaction"


# Records acceleration, displacement, and force time histories along the pier
eval "recorder Node -file $dataDir/pier_accel.out        -time    -nodeRange  [expr int($nodeNum22+$Df/$zSize3)]  [expr int($nodeNum22+$Hpier/$psize+$Df/$zSize3)] -dof 1 2  accel"
eval "recorder Node -file $dataDir/pier_disp.out        -time    -nodeRange   [expr int($nodeNum22+$Df/$zSize3)]  [expr int($nodeNum22+$Hpier/$psize+$Df/$zSize3)] -dof 1 2 3 4 5 6 disp"
eval "recorder Element -file $dataDir/pierforces.out   -time     -eleRange [expr $rec_for_pier+2] [expr int($rec_for_pier+$Hpier/$psize)]    globalForce"

# Records acceleration time history at mid-center of the deck
eval "recorder Node -file $dataDir/Deck_mid_accel.out        -time    -nodeRange  [expr int($nodeNum24-1+($numXele_D/2.0+1)*($numYele_D+1)-$numYele_D)]   [expr int($nodeNum24-1+($numXele_D/2.0+1)*($numYele_D+1))]  -dof 1 2  accel"

# Records displacement time history at the base of the pier 
eval "recorder Node -file $dataDir/pierbase_disp.out        -time    -node  [expr $nodeNum22+2]  -dof 1 2 3 4 5 6 disp"

# Records acceleration time history at the pile cap
eval "recorder Node -file $dataDir/pilecap_accel.out        -time    -node  [expr int($pierbase+2+2.0*$numZnode)]  [expr int($pierbase+2+5.0*$numZnode)] -dof 1 2  accel"

# Records force time history for each 25 piles under the pier --- see my PhD thesis for definition of numbers. 
eval "recorder Element -file $dataDir/pileforces_1.out   -time     -eleRange   $rec_for_pile_1    [expr int($rec_for_pile_1+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_2.out   -time     -eleRange   $rec_for_pile_2    [expr int($rec_for_pile_2+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_3.out   -time     -eleRange   $rec_for_pile_3    [expr int($rec_for_pile_3+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_4.out   -time     -eleRange   $rec_for_pile_4    [expr int($rec_for_pile_4+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_5.out   -time     -eleRange   $rec_for_pile_5    [expr int($rec_for_pile_5+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_6.out   -time     -eleRange   $rec_for_pile_6    [expr int($rec_for_pile_6+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_7.out   -time     -eleRange   $rec_for_pile_7    [expr int($rec_for_pile_7+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_8.out   -time     -eleRange   $rec_for_pile_8    [expr int($rec_for_pile_8+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_9.out   -time     -eleRange   $rec_for_pile_9    [expr int($rec_for_pile_9+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_10.out   -time     -eleRange   $rec_for_pile_10    [expr int($rec_for_pile_10+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_11.out   -time     -eleRange   $rec_for_pile_11    [expr int($rec_for_pile_11+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_12.out   -time     -eleRange   $rec_for_pile_12    [expr int($rec_for_pile_12+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_13.out   -time     -eleRange   $rec_for_pile_13    [expr int($rec_for_pile_13+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_14.out   -time     -eleRange   $rec_for_pile_14    [expr int($rec_for_pile_14+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_15.out   -time     -eleRange   $rec_for_pile_15    [expr int($rec_for_pile_15+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_16.out   -time     -eleRange   $rec_for_pile_16    [expr int($rec_for_pile_16+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_17.out   -time     -eleRange   $rec_for_pile_17    [expr int($rec_for_pile_17+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_18.out   -time     -eleRange   $rec_for_pile_18    [expr int($rec_for_pile_18+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_19.out   -time     -eleRange   $rec_for_pile_19    [expr int($rec_for_pile_19+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_20.out   -time     -eleRange   $rec_for_pile_20    [expr int($rec_for_pile_20+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_21.out   -time     -eleRange   $rec_for_pile_21    [expr int($rec_for_pile_21+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_22.out   -time     -eleRange   $rec_for_pile_22    [expr int($rec_for_pile_22+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_23.out   -time     -eleRange   $rec_for_pile_23    [expr int($rec_for_pile_23+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_24.out   -time     -eleRange   $rec_for_pile_24    [expr int($rec_for_pile_24+$pile_elenum_pier)]     globalForce"
eval "recorder Element -file $dataDir/pileforces_25.out   -time     -eleRange   $rec_for_pile_25    [expr int($rec_for_pile_25+$pile_elenum_pier)]     globalForce"

# Records displacement time history for each 25 piles under the pier --- see my PhD thesis for definition of numbers. 
eval "recorder Node -file $dataDir/piledisp_1.out   -time     -nodeRange   $np_1   [expr int($np_1+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_2.out   -time     -nodeRange   $np_2   [expr int($np_2+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_3.out   -time     -nodeRange   $np_3   [expr int($np_3+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_4.out   -time     -nodeRange   $np_4   [expr int($np_4+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_5.out   -time     -nodeRange   $np_5   [expr int($np_5+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_6.out   -time     -nodeRange   $np_6   [expr int($np_6+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_7.out   -time     -nodeRange   $np_7   [expr int($np_7+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_8.out   -time     -nodeRange   $np_8   [expr int($np_8+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_9.out   -time     -nodeRange   $np_9   [expr int($np_9+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_10.out   -time     -nodeRange   $np_10   [expr int($np_10+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_11.out   -time     -nodeRange   $np_11   [expr int($np_11+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_12.out   -time     -nodeRange   $np_12   [expr int($np_12+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_13.out   -time     -nodeRange   $np_13   [expr int($np_13+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_14.out   -time     -nodeRange   $np_14   [expr int($np_14+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_15.out   -time     -nodeRange   $np_15   [expr int($np_15+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_16.out   -time     -nodeRange   $np_16   [expr int($np_16+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_17.out   -time     -nodeRange   $np_17   [expr int($np_17+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_18.out   -time     -nodeRange   $np_18   [expr int($np_18+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_19.out   -time     -nodeRange   $np_19   [expr int($np_19+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_20.out   -time     -nodeRange   $np_20   [expr int($np_20+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_21.out   -time     -nodeRange   $np_21   [expr int($np_21+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_22.out   -time     -nodeRange   $np_22   [expr int($np_22+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_23.out   -time     -nodeRange   $np_23   [expr int($np_23+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_24.out   -time     -nodeRange   $np_24   [expr int($np_24+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"
eval "recorder Node -file $dataDir/piledisp_25.out   -time     -nodeRange   $np_25   [expr int($np_25+$pile_elenum_pier)] -dof 1 2 3 4 5 6 disp"



#####################################################################################
###### Rigid-Base Boundary - Apply ground motions in the form of displacement  ######
#####################################################################################

# Define ground motions to be applied at the model base ######
set FACT  1.0
timeSeries Path 2000 -fileTime TIME.txt -filePath Disp_X.txt  -factor $FACT
timeSeries Path 3000 -fileTime TIME.txt -filePath Disp_Y.txt  -factor $FACT


pattern Plain 4000 2000 {

for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {
        sp [expr 1+($i-1)*$numYnode*$numZnode+($j-1)*$numZnode] 1 1.0 
      }
    }  
}

pattern Plain 5000 3000 {

for {set i 1} {$i <= $numXnode} {incr i 1} {
  for {set j 1} {$j <= $numYnode} {incr j 1} {

        sp [expr 1+($i-1)*$numYnode*$numZnode+($j-1)*$numZnode] 2 1.0 

      }
    }  
}
########################################################################
#### End of Rigid Base Boundary Condition
########################################################################


########################################################################
#### Compliant Base Boundary Condition - Motion in the form of velocity
########################################################################

# # Create Base Nodes
# source Compliant_Base.tcl

# # User-Defined --- Material properties below model base
# set rockDen 25.0        ;# User-Defined  --- Rock Unit weight
# set rockVs  400.0       ;# User-Defined  --- Rock shear wave velocity
# set rockPoisson 0.25    ;# User-Defined  --- Rock Poisson's ratio

# # Compressive wave (p-wave) --- In case you need to apply motion in vertical direction
# set rockG [expr $rockVs*$rockVs*$rockDen]
# set rockM [expr 2*$rockG*(1-$rockPoisson)/(1-2*$rockPoisson)]
# set rockVp [expr pow($rockM/$rockDen,0.5)]


# # Define Dashpot Material Elements
# set m 0
# for {set i 1} {$i <= $numXnode} {incr i 1} {
  # for {set j 1} {$j <= $numYnode} {incr j 1} {

		# uniaxialMaterial Viscous [expr 500000+$m] [expr $rockDen*$rockVs*$xSize1*$ySize1] 1.0
		# element zeroLength       [expr 600000+$m]  [expr $ID_2+$m]  [expr $ID+$m] -mat [expr 500000+$m]   -dir  1 2
		# set m [expr ($m+1)]
      # }
    # } 
	
# # load pattern for u-rock 
# set FACT  1
# set Gaccel_X 1000
# set Gaccel_Y 2000
# timeSeries Path $Gaccel_X -fileTime Time.txt -filePath Velocity_X.txt  -factor $FACT
# timeSeries Path $Gaccel_Y -fileTime Time.txt -filePath Velocity_Y.txt  -factor $FACT

# set size_average_X 2.50; #User-Defined --- Average element size in direction X
# set size_average_Y 2.50; #User-Defined --- Average element size in direction Y

# pattern Plain 300 1000 {
	# set m 0
		# for {set i 1} {$i <= $numXnode} {incr i 1} {
			# for {set j 1} {$j <= $numYnode} {incr j 1} {
				# load    [expr $ID+$m]    [expr $rockDen*$rockVs*$size_average_X*$size_average_Y] 0.0 0.0 ; # Modify if shaking is applied in other directions.
				# set m [expr ($m+1)]
				# }
			# }

# }

# pattern Plain 301 2000 {
	# set m 0
		# for {set i 1} {$i <= $numXnode} {incr i 1} {
			# for {set j 1} {$j <= $numYnode} {incr j 1} {
				# load    [expr $ID+$m]    0.0 [expr $rockDen*$rockVs*$size_average_X*$size_average_Y] 0.0 ; # Modify if shaking is applied in other directions.
				# set m [expr ($m+1)]
				# }
			# }

# }

########################################################################
#### End of Compliant Base Boundary Condition
########################################################################


##############################################
#### User-Defined Rayleigh Damping
##############################################

set zeta 0.05;		# User-Defined --- Damping Ratio in Decimal
set w1  21.12;		# User-Defined --- Angular Frequency at Frequency 1
set w2  14.96;		# User-Defined --- Angular Frequency at Frequency 2

set a0 [expr $zeta*2.0*$w1*$w2/($w1 + $w2)];	# mass damping coefficient based on first and second modes
set a1 [expr $zeta*2.0/($w1 + $w2)];			# stiffness damping coefficient based on first and second modes

#rayleigh $a0 0.0 $a1 0.0; 

#######################################
###### Perfrom seismic analysis  ######
#######################################

set gamma   0.6
set compdT  0.01  		; #User-Defined --- time-step for transient analysis
set EQduration 30.0		; #User-Defined --- total duration of earthquake

#integrator Newmark  $gamma  [expr pow($gamma+0.5, 2)/4.0]
integrator HHT 0.5 
#numberer RCM
numberer Plain
constraints Penalty 1.0e12 1.0e12
test NormDispIncr 1e-3 50 1
algorithm KrylovNewton
#system SparseGeneral
system SparseSYM
# system Mumps
#system BandGeneral
analysis VariableTransient

 for {set numIncr 1} {$numIncr <= [expr int($EQduration/$compdT)] } {incr numIncr 1} {
        puts "##### Seismic excitation step : $numIncr #####";
        analyze 1 $compdT [expr $compdT/100.0] $compdT  100
        }
set endT [clock seconds]
set ElapsedTime [expr $endT-$startT]
set ElapsedHours [expr int($ElapsedTime/3600.)]
set ElapsedMins [expr int($ElapsedTime/60.-$ElapsedHours*60.)]
set ElapsedSecs [expr int($ElapsedTime-$ElapsedHours*3600-$ElapsedMins*60.)]
puts "Execution time: $ElapsedHours hours and $ElapsedMins minutes and $ElapsedSecs seconds !!"

wipe

