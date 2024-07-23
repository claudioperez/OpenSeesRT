# DispBeamColumnW3D

Warping beam column

[](https://ses.library.usyd.edu.au/bitstream/handle/2123/24045/r917.pdf)

<details><summary>Example script from report</summary>

```tcl
# --------------------------------------------------------------------------------------------------
# Appendix 2. Cantilever beam subject to end torque (Displacement Based Beam-Column Warping Element)
# SET UP ----------------------------------------------------------------------------
wipe; # clear opensees model
source Wsection.tcl; # include definition for I-section
source LibUnits.tcl; # include units
model basic -ndm 3 -ndf 7; # 3 dimensions, 7 dof per node
set dir Cantilever-endtorque
file mkdir $dir; # create data directory
# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1 0 0 0;
node 2 12 0 0 # node#, X Y z
node 3 24 0 0
node 4 36 0 0
node 5 48 0 0
node 6 60 0 0
node 7 72 0 0
node 8 84 0 0
node 9 96 0 0
node 10 108 0 0
node 11 120 0 0
node 12 132 0 0
node 13 144 0 0
node 14 156 0 0
node 15 168 0 0
node 16 180 0 0
node 17 192 0 0
node 18 204 0 0
node 19 216 0 0
node 20 228 0 0
node 21 240 0 0
# Single point constraints -- Boundary Conditions
fix 1 1 1 1 1 1 1 0;
# define material and section
set IDsteel 1; # assign material tag
set Fy 500000; # assign a super large yielding stress
set Es 29000.0;
set Bs 0.000001; # strain-hardening ratio
set R0 15;
set poisson 0.3
set G 11200.0;
set J 5.861
set GJ [expr $G*$J]
set BeamSecTagFiber 1; # assign a tag number to the beam section fiber
set SecTagTorsion 70; # assign a tag number to the torsion information of the beam
set BeamSecTag 3
uniaxialMaterial Steel02 $IDsteel $Fy $Es $Bs; # build steel material
set d [expr 21.62*$in]; # depth
set bf [expr 8.42*$in]; # flange width
set tf [expr 0.93*$in]; # flange thickness
set tw [expr 0.58*$in]; # web thickness
set nfdw 32; # number of fibers along dw
set nftw 4; # number of fibers along tw
set nfbf 32; # number of fibers along bf
set nftf 4; # number of fibers along tf
Wsection $BeamSecTagFiber $IDsteel $d $bf $tf $tw $nfdw $nftw $nfbf $nftf;

#build beam section
uniaxialMaterial Elastic $SecTagTorsion $GJ
section Aggregator $BeamSecTag $SecTagTorsion T -section $BeamSecTagFiber; # add elastic torsion
set numIntgrPts 5; # number of integration points for dis-based element
set BeamTransfTag 1; # associate a tag to beam transformation

# define geometric transformation: performs a corotational geometric
# transformation of beam stiffness and resisting force from the basic system to
# the global-coordinate system
geomTransf Corotational $BeamTransfTag 0 0 1

# Define ELEMENTS -------------------------------------------------------------
element dispBeamColumn 1 1 2 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 2 2 3 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 3 3 4 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 4 4 5 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 5 5 6 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 6 6 7 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 7 7 8 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 8 8 9 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 9 9 10 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 10 10 11 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 11 11 12 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 12 12 13 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 13 13 14 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 14 14 15 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 15 15 16 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 16 16 17 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 17 17 18 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 18 18 19 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 19 19 20 $numIntgrPts $BeamSecTag $BeamTransfTag
element dispBeamColumn 20 20 21 $numIntgrPts $BeamSecTag $BeamTransfTag

# Define RECORDERS -------------------------------------------------------------
# record displacements the free node
recorder Node -file $dir/DFree.out -time -node 21 -dof 1 2 3 4 5 6 7 disp;

# define end torque -------------------------------------------------------------
pattern Plain 1 Linear {
  load 21 0.0 0.0 0.0 1000.0 0.0 0.0 0.0
}

constraints Plain; # how it handles boundary conditions
numberer Plain; # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral; # how to store and solve the system of equations in the analysis
test NormDispIncr 1.0e-8 10 ; # determine if convergence has been achieved at the end of an iteration step
algorithm NewtonLineSearch; # use Newton line search algorithm: updates tangent stiffness at every iteration
integrator DisplacementControl 21 4 0.01; # determine the time step
analysis Static; # define type of analysis static or transient
analyze 100; # perform analysis
puts "Done!"
```

</details>
