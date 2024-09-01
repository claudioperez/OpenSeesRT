import opensees

from opensees import lib, arg
par  = arg.Parameter("Igirz", arg.Num)

con  = lib.ElasticSpring(E=20.0)
lin  = lib.LinearTransform(1)
beam = lib.ElasticBeamColumn3D(nodes=[1, 2], area=20.0, transf=lin, cMass=True)
bea2 = lib.ElasticBeamColumn3D(2, area=20.0, Iz=par, transf=2,   cMass=True, material=con)
lnk  = lib.RigidBeamLink()
abut = lib.ZeroLength3D(1, materials=[0, 2, 1], orientation=[[0.,2],[1,3]])

print(beam.get_cmd())
print(beam(cMass=False,Iz=30).get_cmd())
print(beam(Iy=40).get_cmd())

print(bea2.get_cmd())

print(abut.get_cmd())
#print(abut._args)


assm = opensees.Assembly()

assm.node(1, 0.0, 2.0)
assm.node(2, 0.0, 5.0)
assm.node(3, 5.0, 5.0)

assm.conn(1, "girder",   (1, 2))
assm.conn(2, "abutment", (2, 3))


partials = {
    "girder":   beam,
    "abutment":  abut
}

model = assm.apply(partials)
import emit
model.ndm = 3
model.ndf = 6
print(model.refs)

print(emit.OpenSeesWriter(model).dump())


