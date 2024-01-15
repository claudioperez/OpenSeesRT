#######################################################################################################
#
# procMKPC.tcl
## procedure for evaluating the confined concrete material envelope points based upon the modified
## kent park procedure. The procedure takes in the unconfined concrete and confining steel properties.
## created : NM (nmitra@u.washington.edu) dated: Dec. 2002

#######################################################################################################


proc procMKPC { CUnconfFc CUnconfEc Y Z Cov TSspace TSlength TSFy TSarea Strfactor Lenfactor } {

set CUnconfEcu -0.004;
set SecWid [expr $Lenfactor*$Z]; set SecDep [expr $Lenfactor*$Y]; set cover [expr $Lenfactor*$Cov];
set UFc [expr -$Strfactor*$CUnconfFc]; set Ue0 [expr -$CUnconfEc]; set Uecu [expr -$CUnconfEcu];
set hoopSpc [expr $Lenfactor*$TSspace]; set hoopLngth [expr $Lenfactor*$TSlength];
set hoopFy [expr $Strfactor*$TSFy]; set hoopArea [expr $TSarea*$Lenfactor*$Lenfactor];

# ratio of volume of rectangular steel hoops to volumne of concrete core measured to outside of peripheral hoops
set rhoS [expr ($hoopLngth*$hoopArea)/(($SecWid-2*$cover)*($SecDep-2*$cover)*$hoopSpc)];

# width of concrete core measured to outside of peripheral hoop
set b [expr $SecWid - 2*$cover];
set temp [expr $b/$hoopSpc]
set e50u [expr (3+0.002*$UFc)/($UFc - 1000)]; set e50h [expr 3*$rhoS*pow($temp,0.5)/4];
set Zm [expr 0.5*($UFc-1000)/(3+0.002*$UFc)]; set Z [expr 0.5/($e50u + $e50h - $Ue0)];
set K [expr (1 + $rhoS*$hoopFy/$UFc)];

# unconfined ultimate compressive strength
set UFcu [expr -$UFc*(1-$Zm*($Uecu-$Ue0))/$Strfactor];

#cracking strain in confined concrete
set Ce0 [expr -$K*$Ue0];

# cracking stress in confined concrete
set CFc [expr -$K*$UFc/$Strfactor];

# ultimate stress in confined concrete
set CFcu [expr 0.2*$CFc];

# ultimate strain in confined concrete
set Cecu [expr -(0.8/$Z - $Ce0)];

global concreteProp;
set concreteProp [list $CUnconfFc $CUnconfEc $UFcu $CUnconfEcu $CFc $Ce0 $CFcu $Cecu];

#puts [lindex $concreteProp 0]

return $concreteProp;
}

#===============================================================================================
#===============================================================================================
#===============================================================================================
#===============================================================================================