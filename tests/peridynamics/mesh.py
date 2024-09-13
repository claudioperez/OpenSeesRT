
import meshio
import opensees.openseespy as ops

ndim = 2
totnode = 100
maxfam = 100

model = ops.Model()

model.eval(f"peri init {ndim} {totnode} {maxfam}")

mesh = meshio.read()


for i,node in enumerate(mesh.points):
    x, y, z = ...
    model.eval(f"peri node {i} {x} {y} {z}")


