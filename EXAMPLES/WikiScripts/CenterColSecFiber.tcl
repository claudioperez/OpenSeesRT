# CenterColSecFiber.tcl

# Defines Center Column Section
# Units: kip, in
# KJE, Feb 2003

set h 9.0
set b $h

set beta1 100.0
set beta2 100.0
set spall nospalling
set alphaS 0.015

# Set parameters for fiber section
set nfCore 28
set nfCover   4
set c 1.0
set z [expr $h/2.0-$c]
set y [expr $b/2.0-$c]
set nz [expr -$z]
set ny [expr -$y]
set zc [expr $z+$c]
set yc [expr $y+$c]
set nzc [expr -$zc]
set nyc [expr -$yc]
set fc -3.517
set fc [expr -3.517*$beta2/100.0]
set fcu [expr $fc*$beta1/100.0]

set fy    69.5
set Es    29000.0
set esy   [expr $fy/$Es]
set esu   0.15
set fu [expr $fy+$alphaS*$Es*($esu-$esy)]

	
# Define parameters for elastic section
set Ec 3400
set A  [expr $h*$b]
set Ig [expr $h*$h*$h*$b/12]

set Mcr     228.0
set Kcr     1.2e-4

set My      554.0 ;# for trilinear model (Mmax from UCFyber)
set Ky      6.3e-4 ;# based on UCFyber results

set Ry      [expr $Ky*$L/6.0] ; # yield rotation
set EI      [expr $My/$Ky]
set alpha   0.05

set Ku      0.05
set Mu      [expr $My+$alpha*$EI*($Ku-$Ky)]
set Ru      0.5 ;# this assumes alpha is very small so diff between M-curv alpha and M-rot alpha is negligible

set pinchX  0.5
set pinchY  0.4
set damage1 0.0
set damage2 0.0

set beta    0.4; #only to be used with version 1.3

set Ic      [expr $EI/$Ec]

set slipStiff 91000 ;# see CenterColumnMomCurv.xls
set Rslipy    [expr $My/$slipStiff]
set alphaSlip [expr $alpha/(1.0+(1.0-$alpha)*($slipStiff*$L/6/$EI))]
set Rslipu    [expr $Rslipy+($Mu-$My)/($alphaSlip*$slipStiff)]


# Create axial failure spring
source CenterColAxialSpring.tcl

# Create shear failure spring 
source CenterColShearSpring.tcl

section Aggregator $shearAxialOnlySec $shearTag Vy $axialFailTag P


# Create fiber section
# Define uniaxialMaterials
#                           tag         f'c    epsc   f'cu   epscu
uniaxialMaterial Concrete01  $coreTag   $fc    -0.002  $fcu   -0.0052

if {$spall == "spalling"} {
    uniaxialMaterial Concrete01  $coverTag  $fc    -0.002  0.0    -0.0060
    
} elseif {$spall == "nospalling"} {
    uniaxialMaterial Concrete01  $coverTag  $fc    -0.002  $fcu   -0.0052
    
} else {
    puts stderr "Invalid spalling option: $spall"
}

#                           tag          fy     E      hardening ratio
#uniaxialMaterial Steel02     $steelTag   69.5   29000  $alphaS
uniaxialMaterial Hysteretic $steelTag $fy $esy $fu $esu\
    -$fy -$esy -$fu -$esu 1.0 1.0 0 0

# Define the fiber section
section Fiber $flexSec {
    
    # Define the concrete patch with fibers for unidirectional bending
    patch quadr $coreTag  1 $nfCore $ny $z $ny $nz $y $nz $y $z
    # Define the four cover patches
    patch quadr $coverTag 1 $nfCover $nyc $zc $nyc $nzc $ny $nz $ny $z
    patch quadr $coverTag 1 $nfCover $y $z $y $nz $yc $nzc $yc $zc
    patch quadr $coverTag 1 $nfCore  $ny $nz $nyc $nzc $yc $nzc $y $nz
    patch quadr $coverTag 1 $nfCore  $nyc $zc $ny $z $y $z $yc $zc
    
    
    # Define the reinforcement explicitly using fiber command
    #     yloc  zloc  area matTag
    fiber -3.250  3.250 0.2  $steelTag 0 
    fiber -3.250 -3.250 0.2  $steelTag 0 
    fiber  3.250 -3.250 0.2  $steelTag 0 
    fiber  3.250  3.250 0.2  $steelTag 0 
    fiber -3.187  0.0   0.31 $steelTag 0
    fiber  3.187  0.0   0.31 $steelTag 0
    fiber  0.0   -3.187 0.31 $steelTag 0
    fiber  0.0    3.187 0.31 $steelTag 0 }

# moment-rotation end springs for slip (assume elastic)
uniaxialMaterial Elastic $centerSlipTag $slipStiff

set Acenter $A
set Eccenter $Ec
set Iccenter $Ic

