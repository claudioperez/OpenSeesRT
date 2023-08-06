# Portal frame example that considers three different sections
# (Elastic steel fiber section, Nonlinear steel fiber section, and concrete fiber section) 
# of a force-based beam-column element 

# Vesna Terzic 2011

file mkdir Output; 				   # create data directory
set loadType dynamic;              # pushover or dynamic

foreach rigidConstraint {no yes} {
 
    puts "RIGID CONSTRAINT: $rigidConstraint"
 	
 	foreach secType {concreteFiber steelNonlinearFiber steelElasticFiber} { 
        
		wipe

		# create model builder                                                                                                                                             
		model basic -ndm 2 -ndf 3
		
		set width    360.
		set height   144.

		# create nodes                                                                                                                                                     
		node  1       0.0     0.0
		node  2    $width     0.0
		node  3       0.0 $height
		node  4    $width $height

		# Fix supports at base of columns                                                                                                                                  
		#    tag   DX   DY   RZ                                                                                                                                            
		fix   1     1    1    1
		fix   2     1    1    1
		
		# define mass at nodes 3 and 4
		mass 3 5.0 0.0 0.0
		mass 4 5.0 0.0 0.0
		
		# define constraints
		if {$rigidConstraint == "yes"} {
			equalDOF 3 4 1
		}
		
		# Define materails                                                                                                                       
		# ---------------------------------------------
		set IDcore  1
		set IDcover 2
		set IDsteel 3
		set IDsteelElastic 4
		
		# CONCRETE                     tag      f'c    ec0    f'cu        ecu                                                                                             
		# Core concrete (confined)                                                                                                                                      
		uniaxialMaterial Concrete01  $IDcore   -6.0  -0.004   -3.0     -0.014

		# Cover concrete (unconfined)                                                                                                                                  
		uniaxialMaterial Concrete01  $IDcover  -4.0   -0.002   0.0     -0.005

		# STEEL                                                                                                                                                        
		# Reinforcing steel - nonlinear steel materail                                                                                                                                            
		set fy 60.0;      # Yield stress                                                                                                                               
		set E  29000.0;    # Young's modulus                                                                                                                            
		#                            tag    fy  E0    b  R0  cR1   cR2                                                                                                                       
		uniaxialMaterial Steel02  $IDsteel  $fy $E 0.01  18 0.925  0.15

		# elastic steel materail
		uniaxialMaterial Elastic  $IDsteelElastic  $E
		

		puts "material defined"
		
		# Define cross-section for nonlinear elements                                                                                                                       
		# ---------------------------------------------                                                                                                                       
		
		if { $secType == "concreteFiber" } {
		#FIBER CONCRETE
		source RectangularRCsection2D.tcl

		#                      secID  IDcore  IDcover  IDreinf  d  b cover nfCoreY nfCoreZ nfCoverYlong nfCoverYshort  As   nBarsLayer1 nBarsLayer2 nBarsLayer3
		RectangularRCsection2D   1   $IDcore $IDcover $IDsteel 24. 16.  1.5    20       1         20           2        0.6         3          2           3
		RectangularRCsection2D   2   $IDcore $IDcover $IDsteel 16. 16.  1.5    20       1         20           2        0.6         3          2           3
		}
		
		if { $secType == "steelNonlinearFiber" } {
		#FIBER STEEL
		source WSection.tcl
		#         secID  matID      d      bf   tf    tw   nfdw nftw nfbf nftf
		WSection    1   $IDsteel    23.7 8.97 0.585 0.415    24   1    1    4; # W24x68
		WSection    2   $IDsteel    16.4 16.0 1.89  1.18     24   1    1    4; # W14x257 
		}
		
		if { $secType == "steelElasticFiber" } {
		#FIBER STEEL ELASTIC
		source WSection.tcl
		#         secID  matID              d     bf   tf    tw   nfdw nftw nfbf nftf
		WSection    1   $IDsteelElastic    23.7 8.97 0.585 0.415    24   1    1    4; # W24x68
		WSection    2   $IDsteelElastic    16.4 16.0 1.89  1.18     24   1    1    4; # W14x257
		}
				
		puts "section defined"

		# Define elements                                                                                                                                                  
		# ----------------------                                                                                                                                           

		geomTransf Linear 1
		set np 5
		
		# Create the columns and beam using force-based beam-column elements                                                                                                                   
		#                		tag ndI ndJ nsecs secID transfTag                                                                                                     
		element forceBeamColumn  1   1   3   $np    2       1; # column
		element forceBeamColumn  2   2   4   $np    2       1; # column
		element forceBeamColumn  3   3   4   $np    1       1; # beam

		puts "element defined"
		
		# Define gravity loads                                                                                                                                             
		# --------------------                                                                                                                                             

		# Set a parameter for the axial load                                                                                                                               
		set P 180;                # 10% of axial capacity of columns                                                                                                       

		# Create a Plain load pattern with a Linear TimeSeries
		#                 tag
		timeSeries Linear  1
		#             tag  TStag
		pattern Plain  1     1 {
			# Create nodal loads at nodes 3 & 4                                                                                                                            
			#    nd    FX          FY  MZ                                                                                                                                  
			load  3   0.0  [expr -$P] 0.0
			load  4   0.0  [expr -$P] 0.0
		}

		system BandGeneral
		constraints Transformation
		numberer RCM
		test NormDispIncr 1.0e-12  10 3
		algorithm Newton
		integrator LoadControl 0.1
		analysis Static

		analyze 10

		loadConst -time 0.0
		
		puts "Gravity Done!"
		
		if { $loadType == "pushover" } {
		
			# Define pushover
			pattern Plain 2 1 {
				load  3   1.0  0.0 0.0
				load  4   1.0  0.0 0.0
			}
			
			integrator DisplacementControl 3 1 0.01
			analyze 500; # up to 5 (drift ratio of 3.5%)
			puts "pushover done!"
			
			#Create the output
			set forces [eleResponse 3 forces]
			set axialForce [lindex $forces 0]
			set MzL [lindex $forces 2]
			set MzR [lindex $forces 5]
			set strains [eleResponse 3 basicDeformation]
			set axialStrain [lindex $strains 0]
			set shearForces1 [eleResponse 1 forces]
			set shearForce1 [lindex $shearForces1 0]
			set shearForces2 [eleResponse 2 forces]
			set shearForce2 [lindex $shearForces2 0]
			set disp3 [nodeDisp 3 1]
			set disp4 [nodeDisp 4 1]           
			puts " secType: $secType constraints: $rigidConstraint  axialForce: $axialForce, axialElongation: $axialStrain, shearForce1: $shearForce1, shearForce2: $shearForce2" 
			puts " bendingMoment1: $MzL, bendingMoment2: $MzR"
			puts " dispNode3: $disp3 dispNode4: $disp4"
		
		}
		
		if { $loadType == "dynamic" } {
		
			# Define RECORDERS -------------------------------------------------------------
			recorder EnvelopeElement -file [format "Output/FBeam_%s_%s.out" $secType $rigidConstraint ] -ele 3 force;
			recorder EnvelopeElement -file [format "Output/DBeam_%s_%s.out" $secType $rigidConstraint ] -ele 3 basicDeformation;
			recorder EnvelopeElement -file [format "Output/FCol_%s_%s.out" $secType $rigidConstraint ] -ele 1 2 force;
			
			# DYNAMIC ground-motion analysis -------------------------------------------------------------
			# create load pattern
			if { $secType == "concreteFiber" } {
				set G 386.
			} else {
				set G [expr 3.*386.]
			}
			timeSeries Path 2 -dt 0.005 -filePath A10000.txt -factor $G; # define acceleration vector from file (dt=0.005 is associated with the input file gm)
			pattern UniformExcitation 2 1 -accel 2;		         # define where and how (pattern tag, dof) acceleration is applied
			
			# set damping based on first eigen mode
			set freq [expr [eigen -fullGenLapack  1]**0.5]
			set dampRatio 0.05
			rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]
			
			# create the analysis
			wipeAnalysis;					     # clear previously-define analysis parameters
			constraints Plain;     				 # how it handles boundary conditions
			numberer Plain;					     # renumber dof's to minimize band-width (optimization), if you want to
			system BandGeneral;					 # how to store and solve the system of equations in the analysis
			test NormDispIncr 1.0e-6 10
			algorithm Newton					 # use Linear algorithm for linear analysis
			integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
			analysis Transient;					 # define type of analysis: time-dependent
			analyze 3995 0.01;					 # apply 3995 0.01-sec time steps in analysis
					
			puts "transient Done!"
		
		}

	}
}
wipe