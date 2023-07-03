 #########################################################################################################################
 # Test example for Specimen 3Ca tested by Nisreen Balh (2010) 					  				       #
 # Written: Smail Kechidi, PhD student at University of Blida 1									  	 #
 # mailto: s_kechidi@univ-blida.com															 #
 # PLEASE CITE THIS WITH:																 #
 # Smail Kechidi and Nouredine Bourahla, Deteriorating hysteresis model for cold-formed steel shear wall panel based on  #
 # its physical and mechanical characteristics. Journal of Thin-Walled Structures. 2016; 421-430. 				 #
 # DOI: 10.1016/j.tws.2015.09.022														       #
 # Description: uniaxial material with user defined shear wall panael's physical and mechanical characteristics          #
 # Date: March 12th 2016 																 #
 # Model subjected to The CUREE Cyclic Loading Protocol 											 #
 # File Name: 3Can.tcl 																	 #
 #########################################################################################################################

set N 1.
set sec 1.
set meter 1.
set mm [expr $meter/1000.]
set kN [expr $N*1000.]
set U 1.e2
set E [expr 2.03e+11*$N/pow($meter,2)]
set v 0.3
set G [expr $E/(2*(1+$v))]
set A 10000.0; #rigid 
set Abeam 10000.0; #rigid
set Atruss 10000.0; #rigid

 wipe

 # Create the ModelBuilder object

 model BasicBuilder -ndm 2 -ndf 2

 # Add nodes - command: node nodeId xCrd yCrd
 
    set story_Height 2.44
    set bay_Width 1.22
    node 1 0.0 0.0
    node 2 0.0 $story_Height 
    node 3 $bay_Width $story_Height
    node 4 $bay_Width 0.0
    node 5 [expr $bay_Width/2]  [expr $story_Height/2]
    node 6 [expr $bay_Width/2]  [expr $story_Height/2]

 # Procedure for the loading protocol

 # please keep the follwoing procedures on the same path

 source procRCycDAns.tcl

 # Physical and mechanical characteristics of SWP :

 # Shear Wall Panel's Dimensions :

  set height 2440.0
  set width 1220.0

 # Characteristics and material properties of the steel framing studs :
  set fuf 391
  set fyf 342
  set tf 0.87
  set E 203000.0
  set Af 436.22

 # Characteristics and material properties of sheathing :
  set fus 395.00
  set fys 300.00
  set ts 0.46
  set screw_Spacing 150.0
  set np 1.0

 # Characteristics of the screw fasteners :

  set ds 4.166
  set Vs 1560.0
  set dt 22.2

 # Opening parameters

  set opening_Area 0
  set opening_Length 0

 # material ID

  set matID 1
 
 # add the material to domain through the use of a procedure

 uniaxialMaterial CFSSSWP $matID $height $width $fuf $fyf $tf $Af $fus $fys $ts $np $ds $Vs $screw_Spacing $dt $opening_Area $opening_Length  
 
 uniaxialMaterial Elastic 2 2.03e+11

 element truss 1 1 2 $Atruss 2
 element truss 2 2 3 $Abeam 2
 element truss 3 1 4 $Abeam 2
 element truss 4 3 4 $Atruss 2
 element truss 5 1 5 $Atruss 2
 element truss 6 4 5 $Atruss 2
 element truss 7 2 6 $Atruss 2
 element truss 8 3 6 $Atruss 2

 element zeroLength 9 5 6 -mat 1 -dir 1


 # set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?

 fix 1 1 1
 fix 4 1 1

 pattern Plain 1 Linear {
 load 2 1 0 
 }

 recorder Node -file x3Ca.out -node 3 -dof 1 disp
 recorder Element -file y3Ca.out -ele 9 localForce

 # build the components for the analysis object
 system BandGeneral
 constraints Plain  
 test NormDispIncr 1.0e-8 20
 algorithm Newton
 numberer RCM

 # analysis type used in the procedure is Static

 set peakpts [list   0.0018    0.0018    0.0018    0.0018    0.0018    0.0018    0.0026    0.0020    0.0020    0.0020    0.0020    0.0020    0.0020		0.0035    0.0026    0.0026    0.0026    0.0026    0.0026    0.0026    0.0071    0.0053    0.0053    0.0053    0.0106    0.0081 	0.0081    0.0081    0.0141    0.0106    0.0106    0.0247    0.0187    0.0187    0.0353    0.0265    0.0265]
 set increments 10
 set nodeTag 2
 set dofTag 1

 # start procedure for feeding in Reverse Cyclic loading to the model by Disp. control

 procRCycDAns $increments $nodeTag $dofTag $peakpts
