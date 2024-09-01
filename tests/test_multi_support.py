import opensees as ops
from numpy import sin, linspace, pi

T = 5*pi
t = linspace(0, T, 100)

E = 1.
I = 1.
L = 1.
A = 1.
# k = 60
m = 1.
g = 386.4


beam = ops.element.ElasticBeam2D(E=E, I=I, A=A)

# 3   1   4   2   5
# o---M---o---M---o
# ^       ^       ^
#
data = {
    "nodes": {
        # lumped mass
        1: (  L/2, 0.0),
        2: (3*L/2, 0.0),
        # supports
        3: (  0.0, 0.0),
        4: (    L, 0.0),
        5: (  2*L, 0.0)
    },

    "mass":  {1: m, 2: m},
}

model = ops.model(ndf=3, **data)

model.elem(beam, (3, 1))
model.elem(beam, (1, 4))
model.elem(beam, (4, 2))
model.elem(beam, (2, 5))


pattern = ops.pattern.MultipleSupport([
    (1, 1, ops.pattern.GroundMotion(time=t, accel=(t > T*0.3)*sin(t))),
    (3, 2, ops.pattern.GroundMotion(time=t, accel=(t > T*0.2)*sin(t)))
  ],
)

print(ops.tcl.dumps(pattern))

