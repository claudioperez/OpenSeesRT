import opensees
from numpy import linspace, sin

tag = "1"
fy  = 50e3
E   = 29e6

b   = 0.005

R0  = 18
cR1 = 0.915
cR2 = 0.15

command = f"uniaxialMaterial Steel02 {tag} {fy} {E} {b} {R0} 0.925 0.15"

runtime = opensees.tcl.ModelRuntime(ndm=1, ndf=1)

runtime.eval(command)

strains = 1.2*fy/E*sin(linspace(0, 20, 100))

material = runtime.lift("UniaxialMaterial", tag)

for strain in strains:
    stress = material.getStress(strain, commit=True)
    print(f"{strain}\t{stress}")

