# Eigen analysis of a two-storey one-bay frame; Example 10.5 from "Dynamics of Structures" book by Anil Chopra

# units: kips, in, sec

#       Vesna Terzic, 2010 

#delete all previosly constructed objects
wipe;      

#set input variables
#--------------------

#mass
set m  [expr 100.0/386.0]

#number of modes
set numModes 2

#material
set A 63.41
set I 320.0
set E 29000.0

#geometry
set L 240.
set h  120.

# create data directory
file mkdir modes; 

# define the model
#---------------------------------
#model builder
model BasicBuilder -ndm 2 -ndf 3					

# nodal coordinates:
node 1   0.  0. ;
node 2   $L  0. ;
node 3   0.  $h ;
node 4   $L  $h ;
node 5   0.  [expr 2*$h];
node 6   $L  [expr 2*$h];

# Single point constraints -- Boundary Conditions
fix 1 1 1 1;
fix 2 1 1 1;

# assign mass
mass 3 $m 0. 0. ;
mass 4 $m 0. 0. ;
mass 5 [expr $m/2.] 0. 0. ;
mass 6 [expr $m/2.] 0. 0. ;

# define geometric transformation:
set TransfTag 1;
geomTransf Linear $TransfTag ;

# define elements:
# columns
element elasticBeamColumn 1 1 3 $A $E [expr 2.*$I] $TransfTag;
element elasticBeamColumn 2 3 5 $A $E $I           $TransfTag;
element elasticBeamColumn 3 2 4 $A $E [expr 2.*$I] $TransfTag;
element elasticBeamColumn 4 4 6 $A $E $I           $TransfTag;
# beams
element elasticBeamColumn 5 3 4 $A $E [expr 2.*$I] $TransfTag;
element elasticBeamColumn 6 5 6 $A $E $I           $TransfTag;

# record eigenvectors
#----------------------
for { set k 1 } { $k <= $numModes } { incr k } {
    recorder Node -file [format "modes/mode%i.out" $k] -nodeRange 1 6 -dof 1 2 3  "eigen $k"
}

# perform eigen analysis
#-----------------------------
set lambda [eigen  $numModes];

# calculate frequencies and periods of the structure
#---------------------------------------------------
set omega {}
set f {}
set T {}
set pi 3.141593

foreach lam $lambda {
	lappend omega [expr sqrt($lam)]
	lappend f [expr sqrt($lam)/(2*$pi)]
	lappend T [expr (2*$pi)/sqrt($lam)]
}

puts "periods are $T"

# write the output file cosisting of periods
#--------------------------------------------
set period "modes/Periods.txt"
set Periods [open $period "w"]
foreach t $T {
	puts $Periods " $t"
}
close $Periods

# record the eigenvectors
#------------------------
record 

if 0 {
# create display  for mode shapes
#---------------------------------
#                 $windowTitle $xLoc $yLoc $xPixels $yPixels
# # recorder display "Mode Shape 1"  10    10     500      500     -wipe  
prp $h $h 1;                                         # projection reference point (prp); defines the center of projection (viewer eye)
vup  0  1 0;                                         # view-up vector (vup) 
vpn  0  0 1;                                         # view-plane normal (vpn)     
viewWindow -200 200 -200 200;                        # coordiantes of the window relative to prp  
display -1 5 20;                                     # the 1st arg. is the tag for display mode (ex. -1 is for the first mode shape)
                                                     # the 2nd arg. is magnification factor for nodes, the 3rd arg. is magnif. factor of deformed shape
# # recorder display "Mode Shape 2" 10 510 500 500 -wipe
prp $h $h 1;
vup  0  1 0;
vpn  0  0 1;
viewWindow -200 200 -200 200
display -2 5 20
}

# get values of eigenvectors for translational DOFs
#---------------------------------------------------
set f11 [nodeEigenvector 3 1 1]
set f21 [nodeEigenvector 5 1 1]
set f12 [nodeEigenvector 3 2 1]
set f22 [nodeEigenvector 5 2 1]
puts "eigenvector 1: [list [expr {$f11/$f21}] [expr {$f21/$f21}] ]"
puts "eigenvector 2: [list [expr {$f12/$f22}] [expr {$f22/$f22}] ]"














