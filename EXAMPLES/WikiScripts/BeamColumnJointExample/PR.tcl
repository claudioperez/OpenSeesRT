
############################################################################
# Test example for BEAM COLUMN ELEMENT JOINT ------- PARK RUITONG TEST SPECIMEN Unit 1
# Written: N.Mitra
# Description: 4 noded 12 dof element having 12 springs and a shear panel
# Date: Feb 16 2003
## Model consisting of a crucifix with beams and columns and a joint
## File Name: PR1.tcl
# refer to Beam-Column-Joint Element.doc for full explanation of the parameters

################################################################################

#create the ModelBuilder object

model BasicBuilder -ndm 2 -ndf 3

## unit name ---- PR
set fName "PR";

source ProcMKPC2.tcl
puts "a"
#tsource procUniaxialPinching.tcl
source ProcRC2.tcl

puts "a"
# all dimensions are in here as MPa (conversion factors are required in certain places)
set Strfactor 145; set Lenfactor [expr 1/25.4];

## Y taken as the inplane dim. against which the bending takes place
set colY 406; set colZ 305;
set bmY 457; set bmZ 229;

# covers
set colCov 43; set bmCov1 42; set bmCov2 33; set bmCov $bmCov1;

# y,z,x dimension of the joint respectively
set JointWidth [expr $colY]; set JointHeight [expr $bmY]; set JointDepth $colZ ;
set BeamLengthIn 645; set BeamLengthOut 1271; set ColumnLengthClear 1008;
set JointVolume [expr $JointWidth*$JointHeight*$JointDepth];

############################################################
################# material properties of column section ############################################################

set CUnconfFc -45.9; set CUnconfEc -0.002;
set CTSspace 60; set CTSlength 1853.53; set CTSFy 282; set CTSarea 28.3;
set CFy 498.0; set CEs 196600.0; set CsHratio 0.004216; set CAs 201.06;

procMKPC $CUnconfFc $CUnconfEc $colY $colZ $colCov $CTSspace $CTSlength $CTSFy $CTSarea $Strfactor $Lenfactor

set CUnconfFcu [lindex $concreteProp 2]; set CUnconfEcu [lindex $concreteProp 3];
set CConfFc [lindex $concreteProp 4]; set CConfEc [lindex $concreteProp 5];
set CConfFcu [lindex $concreteProp 6]; set CConfEcu [lindex $concreteProp 7];

#########################################################
########################### material properties of beam section #####################################################

set BUnconfFc -45.9; set BUnconfEc -0.002;
set BTSspace 80; set BTSlength 1036; set BTSFy 282; set BTSarea 28.3;
set BFy 294.0; set BEs 210400.0; set BAs 201.06; set BsHratio 0.002322;

procMKPC $BUnconfFc $BUnconfEc $bmY $bmZ $bmCov $BTSspace $BTSlength $BTSFy $BTSarea $Strfactor $Lenfactor

set BUnconfFcu [lindex $concreteProp 2]; set BUnconfEcu [lindex $concreteProp 3];
set BConfFc [lindex $concreteProp 4]; set BConfEc [lindex $concreteProp 5];
set BConfFcu [lindex $concreteProp 6]; set BConfEcu [lindex $concreteProp 7];

####################################################################
######################### details for the material models of bar slip of the beam ####################################

set bs_fc [expr -$BUnconfFc]; set bs_fs $BFy; set bs_es $BEs; set bs_fsu 434; set bs_dbar 16; set bs_esh [expr $BsHratio*$BEs];
set bs_wid $colZ; set bs_dep $bmY;
set bsT_nbars 5; set bsB_nbars 2;
set bs_ljoint $colY;

################################################################################
########################## details for the material models of bar slip of the column ##################################

set cs_fc [expr -$CUnconfFc]; set cs_fs $CFy; set cs_es $CEs; set cs_fsu 660; set cs_dbar 16; set cs_esh [expr $CsHratio*$CEs];
set cs_wid $colZ; set cs_dep $colY;
set cs_nbars 3;
set cs_ljoint $bmY;

#####################################################################
############### add nodes - command: node nodeId xCrd yCrd #######################################################

node 1 0.0 0.0
node 2 0.0 $ColumnLengthClear
node 3 [expr -$BeamLengthOut-$BeamLengthIn-$JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 4 [expr -$BeamLengthIn-$JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 5 [expr -$JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 6 [expr $JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 7 [expr $BeamLengthIn+$JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 8 [expr $BeamLengthOut+$BeamLengthIn+$JointWidth/2] [expr $ColumnLengthClear+$JointHeight/2]
node 9 0.0 [expr $ColumnLengthClear+$JointHeight]
node 10 0.0 [expr 2*$ColumnLengthClear+$JointHeight]

# add material Properties - command: uniaxialMaterial matType matTag ...
#command: uniaxialMaterial Elastic tag? E?

uniaxialMaterial Elastic 1 10000000000.0

#####################################################################
############## inelastic beam column elements
##########################################################

uniaxialMaterial Concrete01 10 $BUnconfFc $BUnconfEc $BUnconfFcu $BUnconfEcu
uniaxialMaterial Concrete01 20 $BConfFc $BConfEc $BConfFcu $BConfEcu
uniaxialMaterial Steel02 30 $BFy $BEs $BsHratio 18.5 0.925 0.15 0.0 0.4 0.0 0.5
uniaxialMaterial Concrete01 40 $CUnconfFc $CUnconfEc $CUnconfFcu $CUnconfEcu
uniaxialMaterial Concrete01 50 $CConfFc $CConfEc $CConfFcu $CConfEcu
uniaxialMaterial Steel02 60 $CFy $CEs $CsHratio 18.5 0.925 0.15 0.0 0.4 0.0 0.5

########### for columns ///////////////////////////////////////////////////////////////////////

set z [expr $colZ/2.0]; set y [expr $colY/2.0];

section Fiber 1 {
patch rect 50 8 1 [expr $colCov-$y] [expr $colCov-$z] [expr $y-$colCov] [expr $z-$colCov]
patch rect 40 2 1 [expr -$y] [expr $colCov-$z] [expr $colCov-$y] [expr $z-$colCov]
patch rect 40 2 1 [expr $y-$colCov] [expr $colCov-$z] [expr $y] [expr $z-$colCov]
patch rect 40 8 1 [expr -$y] [expr -$z] [expr $y] [expr $colCov-$z]
patch rect 40 8 1 [expr -$y] [expr $z-$colCov] [expr $y] [expr $z]

layer straight 60 3 $CAs [expr $y-$colCov] [expr $colCov-$z] [expr $y-$colCov] [expr $z-$colCov]
layer straight 60 2 $CAs 0.0 [expr $colCov-$z] 0.0 [expr $z-$colCov]
layer straight 60 3 $CAs [expr $colCov-$y] [expr $colCov-$z] [expr $colCov-$y] [expr $z-$colCov]
}

#################### for beams ///////////////////////////////////////////////////////////////////

set z [expr $bmZ/2.0]; set y [expr $bmY/2.0];

section Fiber 2 {
patch rect 20 8 1 [expr $bmCov1-$y] [expr $bmCov1-$z] [expr $y-$bmCov1] [expr $z-$bmCov1]
patch rect 10 2 1 [expr -$y] [expr $bmCov1-$z] [expr $bmCov1-$y] [expr $z-$bmCov1]
patch rect 10 2 1 [expr $y-$bmCov1] [expr $bmCov1-$z] [expr $y] [expr $z-$bmCov1]
patch rect 10 8 1 [expr -$y] [expr -$z] [expr $y] [expr $bmCov1-$z]
patch rect 10 8 1 [expr -$y] [expr $z-$bmCov1] [expr $y] [expr $z]

layer straight 30 3 $BAs [expr $y-$bmCov1] [expr $bmCov1-$z] [expr $y-$bmCov1] [expr $z - $bmCov1]
layer straight 30 2 $BAs [expr $y-$bmCov1-$bmCov2] [expr $bmCov1-$z] [expr $y-$bmCov1-$bmCov2] [expr $z - $bmCov1]
layer straight 30 2 $BAs [expr $bmCov1-$y] [expr $bmCov1-$z] [expr $bmCov1-$y] [expr $z - $bmCov1]
}

##############/////////////////////////////////////////////////////////////////////////////////////

## add geometric transformation -command: geomTransf transfType ...

## geomTransf Linear tag?

geomTransf Linear 1
geomTransf Linear 2

element nonlinearBeamColumn 1 1 2 5 1 2
element nonlinearBeamColumn 2 9 10 5 1 2
element nonlinearBeamColumn 3 3 4 3 2 1
element nonlinearBeamColumn 4 4 5 2 2 1
element nonlinearBeamColumn 5 6 7 2 2 1
element nonlinearBeamColumn 6 7 8 3 2 1

##### end element formation as well as material defination for beams and columns ######################

##########################################################################################

# for beam bottom

set matID1 21
set matID2 22

uniaxialMaterial BarSlip $matID1 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsB_nbars $bs_wid $bs_dep strong beamBot
uniaxialMaterial BarSlip $matID2 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsB_nbars $bs_wid $bs_dep strong beamBot

## %%%%%%%%%%%%%%% equivalent statement can be made in other way
#uniaxialMaterial BarSlip $matID1 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsB_nbars $bs_wid $bs_dep 1.0 strong beamBot damage
#uniaxialMaterial BarSlip $matID2 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsB_nbars $bs_wid $bs_dep 1.0 strong beamBot damage

# for beam top
set matID3 31
set matID4 32

uniaxialMaterial BarSlip $matID3 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsT_nbars $bs_wid $bs_dep strong beamTop
uniaxialMaterial BarSlip $matID4 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsT_nbars $bs_wid $bs_dep strong beamTop

## %%%%%%%%%%%%%%% equivalent statement can be made in other way

#uniaxialMaterial BarSlip $matID3 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsT_nbars $bs_wid $bs_dep 1.0 strong beamTop damage
#uniaxialMaterial BarSlip $matID4 $bs_fc $bs_fs $bs_es $bs_fsu $bs_esh $bs_dbar $bs_ljoint $bsT_nbars $bs_wid $bs_dep 1.0 strong beamTop damage

# for columns
set matID5 41
set matID6 42
set matID7 43
set matID8 44

uniaxialMaterial BarSlip $matID5 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep strong column
uniaxialMaterial BarSlip $matID6 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep strong column
uniaxialMaterial BarSlip $matID7 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep strong column
uniaxialMaterial BarSlip $matID8 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep strong column

## %%%%%%%%%%%%%%% equivalent statement can be made in other way

#uniaxialMaterial BarSlip $matID5 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep 1.0 strong column damage
#uniaxialMaterial BarSlip $matID6 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep 1.0 strong column damage
#uniaxialMaterial BarSlip $matID7 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep 1.0 strong column damage
#uniaxialMaterial BarSlip $matID8 $cs_fc $cs_fs $cs_es $cs_fsu $cs_esh $cs_dbar $cs_ljoint $cs_nbars $cs_wid $cs_dep 1.0 strong column damage

#################### end material formation for bar slip ########################################################

############################ material for shear panel #############################################################
## Positive/Negative envelope Stress
set p1 2.1932; set p2 4.0872; set p3 4.4862; set p4 [expr $p3*1e-3];

## stress1 stress2 stress3 stress4
set pEnvStrsp [list [expr $p1*$JointVolume] [expr $p2*$JointVolume] [expr $p3*$JointVolume] [expr $p4*$JointVolume]]
set nEnvStrsp [list [expr -$p1*$JointVolume] [expr -$p2*$JointVolume] [expr -$p3*$JointVolume] [expr -$p4*$JointVolume]]

## Positive/Negative envelope Strain
## strain1 strain2 strain3 strain4

set pEnvStnsp [list 0.0002 0.004465 0.0131 0.0269]
set nEnvStnsp [list -0.0002 -0.004465 -0.0131 -0.0269]

## Ratio of maximum deformation at which reloading begins
## Pos_env. Neg_env.
set rDispsp [list 0.25 0.25]

## Ratio of envelope force (corresponding to maximum deformation) at which reloading begins

### Pos_env. Neg_env.
set rForcesp [list 0.15 0.15]


## Ratio of monotonic strength developed upon unloading
### Pos_env. Neg_env.

set uForcesp [list 0.0 0.0]


## Coefficients for Unloading Stiffness degradation

## gammaK1 gammaK2 gammaK3 gammaK4 gammaKLimit

set gammaKsp [list 1.13364492409642 0.0 0.10111033064469 0.0 0.91652498468618]

#set gammaKsp [list 0.0 0.0 0.0 0.0 0.0]

#### Coefficients for Reloading Stiffness degradation
### gammaD1 gammaD2 gammaD3 gammaD4 gammaDLimit

set gammaDsp [list 0.12 0.0 0.23 0.0 0.95]

#set gammaDsp [list 0.0 0.0 0.0 0.0 0.0]
#### Coefficients for Strength degradation
### gammaF1 gammaF2 gammaF3 gammaF4 gammaFLimit

set gammaFsp [list 1.11 0.0 0.319 0.0 0.125]
#set gammaFsp [list 0.0 0.0 0.0 0.0 0.0]

set gammaEsp 10.0

uniaxialMaterial Pinching4 5 [lindex $pEnvStrsp 0] [lindex $pEnvStnsp 0] \
[lindex $pEnvStrsp 1] [lindex $pEnvStnsp 1] [lindex $pEnvStrsp 2] \
[lindex $pEnvStnsp 2] [lindex $pEnvStrsp 3] [lindex $pEnvStnsp 3] \
[lindex $nEnvStrsp 0] [lindex $nEnvStnsp 0] \
[lindex $nEnvStrsp 1] [lindex $nEnvStnsp 1] [lindex $nEnvStrsp 2] \
[lindex $nEnvStnsp 2] [lindex $nEnvStrsp 3] [lindex $nEnvStnsp 3] \
[lindex $rDispsp 0] [lindex $rForcesp 0] [lindex $uForcesp 0] \
[lindex $rDispsp 1] [lindex $rForcesp 1] [lindex $uForcesp 1] \
[lindex $gammaKsp 0] [lindex $gammaKsp 1] [lindex $gammaKsp 2] [lindex $gammaKsp 3] [lindex $gammaKsp 4] \
[lindex $gammaDsp 0] [lindex $gammaDsp 1] [lindex $gammaDsp 2] [lindex $gammaDsp 3] [lindex $gammaDsp 4] \
[lindex $gammaFsp 0] [lindex $gammaFsp 1] [lindex $gammaFsp 2] [lindex $gammaFsp 3] [lindex $gammaFsp 4] \
$gammaEsp energy

#################### end material formation for shear panel ###########################################

##element BeamColumnJoint tag? iNode? jNode? kNode? lNode? matTag1? matTag2? matTag3? matTag4?
## matTag5? matTag6? matTag7? matTag8? matTag9? matTag10? matTag11? matTag12? matTag13?
## <element Height factor?> <element Width factor?>
## please note: the four nodes are in anticlockwise direction around the element
## requires material tags for all 13 different components within the element.
## the first 12 being that of spring and the last of the shear panel

element beamColumnJoint 7 2 6 9 5 41 42 1 21 31 1 43 44 1 22 32 1 5

#element beamColumnJoint 7 2 6 9 5 1 1 1 1 1 1 1 1 1 1 1 1 1

## %%%%%%%%%%%%%%% equivalent statement can be made in other way

#element beamColumnJoint 7 2 6 9 5 41 42 1 21 31 1 43 44 1 22 32 1 5 1.0 1.0
#element beamColumnJoint 7 2 6 9 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1.0 1.0

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
fix 1 1 1 0
fix 2 0 0 0
fix 3 0 1 0
fix 4 0 0 0
fix 5 0 0 0
fix 6 0 0 0
fix 7 0 0 0
fix 8 0 1 0
fix 9 0 0 0
fix 10 0 0 0

pattern Plain 2 Linear {
load 4 0 -55000 0 -const
load 7 0 -55000 0 -const
}

system ProfileSPD
constraints Plain
integrator LoadControl 0 1 0 0
test NormDispIncr 1e-8 150
algorithm Newton
numberer RCM
analysis Static
analyze 1
loadConst -time 0.0

pattern Plain 1 Linear {
#load nd? Fx? Fy? Mz?
load 10 1 0 0
}

set rbbt "_RBbt"; set rbtp "_RBtp"; set dlcbr "_DLCbr"; set sp "_Sp"; set jdf "_Jdf";
set lbbt "_LBbt"; set lbtp "_LBtp"; set drcbr "_DRCbr"; set ulcbr "_ULCbr"; set urcbr "_URCbr";
set RBbt [concat $fName$rbbt]; set RBtp [concat $fName$rbtp]; set Sp [concat $fName$sp];
set LBbt [concat $fName$lbbt]; set LBtp [concat $fName$lbtp]; set DLCbr [concat $fName$dlcbr];
set DRCbr [concat $fName$drcbr]; set URCbr [concat $fName$urcbr]; set ULCbr [concat $fName$ulcbr];
set Jdf [concat $fName$jdf];


recorder Node -file $fName.out -time -node 10 -dof 1 disp;
#recorder Node $fName.out disp -load -node 10 -dof 1
#recorder Element 7 -file $RBbt.out node2BarSlipB stressStrain
#recorder Element 7 -file $RBtp.out node2BarSlipT stressStrain
#recorder Element 7 -file $LBbt.out node4BarSlipB stressStrain
#recorder Element 7 -file $LBtp.out node4BarSlipT stressStrain
#recorder Element 7 -file $DLCbr.out node1BarSlipL stressStrain
#recorder Element 7 -file $DRCbr.out node1BarSlipR stressStrain
#recorder Element 7 -file $ULCbr.out node3BarSlipL stressStrain
#recorder Element 7 -file $URCbr.out node3BarSlipR stressStrain
#recorder Element 7 -file $Sp.out shearpanel stressStrain
#recorder Element 7 -file $Jdf.out deformation

set peakpts [list 0.1 10 10 30 30 45 45 60 60 75 75 90 90 105 105]
set increment 1
set nodeTag 10
set dofTag 1

procRC $increment $nodeTag $dofTag $peakpts

# print the results at node and at all elements

print node

#print element




set disp [nodeDisp 10 1]

puts "disp = $disp"
