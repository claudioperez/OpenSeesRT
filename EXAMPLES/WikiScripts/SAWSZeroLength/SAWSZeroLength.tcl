## Units kips-inches
##SAWS Model test file - static
##----------------------------------------------------------
model      BasicBuilder   -ndm 2        -ndf 3
## Nodes
node 1 0.0   0.0
node 2 0.0   0.0
## Mass
mass 2 1.0   0.0    0.0
#Boundary Conditions
fix 1 1 1 1
fix 2 0 1 1
## Materials  
set F0 [expr 8.0*0.225]  ;# kip
set FI [expr 1.20*0.225]  ;# kip
set DU [expr 15.0/25.4]   ;# in
set S0 [expr 5.0*0.225*25.4]  ;# kip/in
set R1 0.058
set R2 -0.050
set R3 1.00
set R4 0.020
set alph 0.60
set bet  1.10 

#set F0 5.28
#set FI 3.28
#set DU 0.84
#set S0 87.84
#set R1 0.12
#set R2 -0.02
#set R3 1.00
#set R4 0.08
#set alph 0.7
#set bet  1.10 

#uniaxialMaterial SAWS 1 $F0 $FI $DU $S0 $R1 $R2 $R3 $R4 $alph $bet
uniaxialMaterial SAWS  1 15.7998483699773 0.545094768764215 1.04095 159.837064220183 0.102201841215193 -0.0324118361176701 1 0.0692552595908408 0.8 1.1
uniaxialMaterial Elastic 2 1e10
## Transformation
geomTransf Linear 1 
## Define Model
element zeroLength 1 1 2 -mat 1 2   -dir 1 2 


