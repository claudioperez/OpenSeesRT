    #########################################################
    # Cyclic analysis of a Lehman's Column 415 (PEER 1998/01)
    #########################################################

    # Written: Vesna Terzic (vesna@berkeley.edu)
    # Created: 12/2011

    #units: kip,in

    wipe

    source LibUnits.tcl

    set eleType 1;  # force-based = 1
                    # displacement-based = 2

    #number of finite elements to model the column and the number of integration points per element
    if { $eleType == 1 } {
        set NoEle 1
        set nIP 5
    } elseif { $eleType == 2 } {
        set NoEle 4
        set nIP 3
    }

    # ------------------------------
    # Start of model generation
    # ------------------------------
    # create ModelBuilder (with 2-dimensions and 3 DOF/node)
     model BasicBuilder -ndm 2 -ndf 3


    #input parameters
    set HCol [expr 96.*$in]; #column height
    set DCol [expr 24.*$in]; #column diameter
    set clearCover [expr 0.75*$in]; #clear cover of concrete

    #derived quantities
    set HEle [expr $HCol/$NoEle]; #height of one finite element
    set ACol [expr $DCol**2*$PI/4.0]; #area of the column cross-section

    # Define geometry for model
    # -------------------------
    #node $tag $xCrd $yCrd
     node  1     0.0  0.0
    for { set i 0 } { $i < $NoEle } { incr i } {
        node  [expr $i+2]   0.0  [expr ($i+1)*$HEle]
    }

    # set the boundary conditions
    # fix tag DX DY Rot
     fix  1   1  1  1

    # Define coordinate transformation
    # --------------------------------
    set transfTag 1
    # geomTransf Linear $transfTag
    # geomTransf PDelta $transfTag
    geomTransf Corotational $transfTag

    #------------------------------------------------------
    # Define uniaxial materials
    #(the three materails will be defined: one for unconfined concrete cover, one for concrete confined core, and one for reinforcement)
    #--------------------------------------------------------

    # define material tags
    set coreTag 1
    set coverTag 2
    set steelTag 3

    #define longitudinal reinforcement
    set barArea [expr 0.31*$in2]; #area of longitudinal bar #5
    set db [expr 0.625*$in]; #longitudinal bar diameter
    set numBars 22; #number of longitudinal bars
    set fy [expr 70.*$ksi]; #yield strength of longitudinal bars
    set Es [expr 29000.*$ksi]; #modulus of elasticity of steel
    set Esf [expr $Es*0.025]; # tangent at initial strain hardening (calibrated from the coupon test)


    #define transverse reinforcement
    set dh [expr 0.25*$in]; #diameter of the spiral #2
    set NoHoops 1; #number of hoops in the bundle
    set Asp1 [expr 0.0491*$in2]; #area of transverse reinforcement bar
    set stran [expr 1.25*$in]; #the centerline distance between spirals along the height of the column
    set fyh [expr 96.6*$ksi]; #yield strenght of the hoop

    #------------------------------------------
    #  define steel model
    #------------------------------------------

    uniaxialMaterial Steel02   [expr $steelTag+1]  $fy $Es 0.025 18.0 0.925 0.15; # Menegotto-Pinto uniaxial steel model (coefficients by Terzic, 2010)
    uniaxialMaterial MinMax     $steelTag  [expr $steelTag+1]  -min -0.080 -max 0.080

    #-------------------------------
    #define concrete model
    #-------------------------------

    # Plain (unconfined concrete)

    #compression
    set fc [expr 4.4*$ksi]; # compressive strenght of plain concrete on the day of the test (data from the concrete cylinder tests)
    set eps0 0.002; #strain that corresponds to fc' (Caltrans SDC)
    set epss 0.005; #ultimate strain for unconfined concrete (Caltrans SDC)
    set Ec [expr 57000.*sqrt($fc*1000.)/1000.]; # the formula from ACI building code

    #Confined concrete

    #compression
    # Mander's equations for calculating confined concrete compressive strength and corresponding strain (Mander, 1988)
    set sprime [expr $stran-$NoHoops*$dh]; #clear distance between spirals
    set ds [expr $DCol-2.0*$clearCover-$dh]; #diameter of spirals between spiral bar centers (diameter of the confined core)
    set Asp [expr $Asp1*$NoHoops]; #total area of transverse reinforcement in the bundle

    set As [expr $barArea*$numBars]; #total area of the longitudinal steel in the section
    set Ac [expr $ds**2*$PI/4.0]; #area of core of section
    set rho_cc [expr $As/$Ac]; #ratio of area of longitudinal reinf. to area of core of section

    set ke [expr (1.0 - $sprime/2.0/$ds)/(1.0-$rho_cc)]; #confinement effectiveness coeficient = Ae/Acc
    set rho_t [expr 4.0*$Asp/$ds/$stran]; # ratio of transverse reinforcement
    set fl [expr 1.0/2.0*$ke*$rho_t*$fyh]; #effective lateral confining stress on the concrete

    set fcc [expr $fc*(-1.254+2.254*sqrt(1.0+7.94*$fl/$fc)-2.0*$fl/$fc)]; #confined concrete compressive strength
    set epscc [expr $eps0*(1.0+5.0*($fcc/$fc-1.0))]; #strain that corresponds to fcc'

    # ultimate stress and strain
    set ecr [expr $Ec/($Ec-$fcc/$epscc)]; #r factor (Eq. 6 in Mander, 1988)
    set epscu [expr 0.004+0.14*$fyh/$fc*$rho_t]; #ultimate strain (by Dawn Lehman,1998, PEER 1998/01)
    set fcu [expr $fcc*$epscu/$epscc*$ecr/($ecr-1.0+pow($epscu/$epscc,$ecr))]; #strength that corresponds to ultimate strain (Mander, 1988)

    uniaxialMaterial Concrete02 $coverTag  -$fc   [expr -2.0*$fc/$Ec]    0.0  -$epss 0.1 [expr -0.04*$fc] [expr 0.04*$fc/(2.0*$fc/$Ec)] ; #  plain concrete
    uniaxialMaterial Concrete02 $coreTag   -$fcc  [expr -2.0*$fcc/$Ec] -$fcu  -$epscu 0.1 [expr -0.04*$fcc] [expr 0.04*$fcc/(2.0*$fcc/$Ec)]; # confined concrete

    #-----------------------
    #Define Section
    #-----------------------
    set secnTag 1

    set nfCoreT 12
    set nfCoreR 8
    set nfCoverT 12
    set nfCoverR 2
    set total_conc_fibers [expr $nfCoreT*$nfCoreR+$nfCoverT*$nfCoverR]

    set ro [expr $DCol/2.0]; #radius of the column cross-section
    set rl [expr $ro-$clearCover-$dh-($db/2.0)]; #distance from the column centroid to the centroid of the long. bar
    set ri [expr ($ds-$sprime/4.)/2.]; #radius of the effectively confined core + s'/8

    # define circular fiber section
    section fiberSec $secnTag -GJ 1e8 {
        # Define the core patch
        patch circ $coreTag $nfCoreT $nfCoreR 0 0 0 $ri 0.0 360.0
        # Define the cover patch
        patch circ $coverTag $nfCoverT $nfCoverR 0 0 $ri $ro 0.0 360.0
        # Define the reinforcing layer
        set theta [expr 360.0/$numBars]
        layer circ $steelTag $numBars $barArea 0 0 $rl [expr $theta/2.] [expr 360.0-$theta/2.]
    }

    # Caltrans column shear capacity in English units (Chapter 3 of Caltrans SDC)
    set poisson 0.20
    set G [expr $Ec/2.0/(1.0+$poisson)]

    set fv1 [expr $rho_t*$fyh/0.15+3.67-3.0]
    if { $fv1 < 0.3 } {
        set fv1 0.3
    }
    if { $fv1 > 3 } {
        set fv1 3
    }
    set ALR [expr 0.0475*$ACol*$fc*1000./2000./$ACol]
    set fv2 [expr 1.0+$ALR]
    if { $fv2 > 1.5 } {
        set fv2 1.5
    }
    set vc [expr $fv1*$fv2*sqrt($fc*1000.0)]
    if { $vc > [expr 4.*sqrt($fc*1000.0)] } {
        set vc [expr 4.*sqrt($fc*1000.0)]
    }
    set Vc [expr 0.8*$ACol*$vc/1000.0]

    set Vs [expr $Asp1*$PI/2.0*$fyh*$ds/$stran]
    set Vn [expr $Vc+$Vs]
    set Vst [expr 3.0/4.0*$G*$ACol]
    set gam_y [expr $Vn/$Vst]

    #define shear force-deformation relationship
    #                                matTag         Fy     E0     b
    uniaxialMaterial Steel01  [expr $steelTag+6]    $Vn   $Vst  1.0e-3

    # Aggregate shear to the RC section
    set secTag [expr $secnTag+1]
    section Aggregator $secTag  [expr $steelTag+6]  Vy   -section $secnTag

    #------------------------------
    # Define elements
    # -----------------------------

    for { set i 0 } { $i < $NoEle } { incr i } {
        if { $eleType == 1 } {
            element forceBeamColumn  [expr $i+1]   [expr $i+1] [expr $i+2]  $nIP $secTag $transfTag
        } elseif { $eleType == 2 } {
            element dispBeamColumn  [expr $i+1]   [expr $i+1] [expr $i+2]  $nIP $secTag $transfTag
        }
    }

    #---------------------------------------------
    # Define Gravity Load
    #---------------------------------------------
    set IDctrlNode [expr $NoEle+1]

    pattern Plain 1 Linear {
       load $IDctrlNode 0. -147.0 0.;    # node#, FX FY MZ
    }
    constraints Plain;                     # how it handles boundary conditions
    numberer Plain;                        # renumber dof's to minimize band-width (optimization), if you want to
    system BandGeneral;                    # how to store and solve the system of equations in the analysis
    test NormDispIncr 1.0e-8 6 ;         # determine if convergence has been achieved at the end of an iteration step
    algorithm Newton;                    # use Newton's solution algorithm: updates tangent stiffness at every iteration
    integrator LoadControl 0.1;            # determine the next time step for an analysis, # apply gravity in 10 steps
    analysis Static                        # define type of analysis
    analyze 10;                            # perform gravity analysis
    loadConst -time 0.0;                # hold gravity constant and restart time

    #--------------------------------
    # Cyclic load pattern
    #--------------------------------

    # write displacement history into a file
    set fileu [open out/displacement.txt "w"]
    close $fileu

    set ductility 7
    set uy 1.0
    set n 36
    set dhpret [list 0.06 0.06 0.06 0.15 0.15 0.15 0.3 0.3 0.3 0.75 0.75 0.75 1.00 1.00 1.00]
    #set dhpost [lappend dhpret 1.5 1.5 1.5 0.5 2.0 2.0 2.0 0.65 3.0 3.0 3.0 1.0 5 5 5 1.66 7 7 7]
    set dhpost [lappend dhpret 1.5 1.5 1.5 0.5 2.0 2.0 2.0 0.65 3.0 3.0 3.0 1.0 5 5 5]
    set dhtot $dhpost
    set running 0.0

    set cycles [llength $dhtot]
    source SingleCycle.tcl
    for { set k 0 } { $k < $cycles } { incr k } {
        set cycmax [expr $uy*[lindex $dhtot $k]]
        set thist [singlecycle $cycmax $n $running]
        set running [expr $running + $thist]
    }

    set anpts [expr (4*$n-3)*$cycles]
    set dt [expr $running/($anpts-1)]
    puts "pts=$anpts, dt=$dt, tfinal=$running"

    # define time series
    timeSeries Path 2 -dt $dt -filePath out/displacement.txt
    set IDctrlDOF 1;

    # define load pattern
    pattern Plain 2 2 {
        sp $IDctrlNode $IDctrlDOF 1.0
    }

    # define recorders for displacement and force
    recorder Node -file out/Disp.out -time -node [expr $NoEle +1] -dof 1 disp; #records displacement at the top node (TN)
    recorder Node -file out/Force.out -time -node 1 -dof 1 reaction; #records displacement at the top node (TN)

    # cyclic analysis objects
    constraints Penalty 1.0e14 1.0e14
    integrator LoadControl  $dt
    numberer Plain
    system BandGeneral
    test NormDispIncr      1.0e-7    10        0
    algorithm Newton
    analysis Static

    analyze $anpts




