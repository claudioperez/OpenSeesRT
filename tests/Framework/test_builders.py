import opensees.tcl as pyg3

a, b = [pyg3.SafeBuilder(ndm=2, ndf=3)  for i in [1,2]]

c, d = [pyg3.BasicBuilder(ndm=2, ndf=3) for i in [1,2]]


assert a._domain != b._domain

assert c._domain == d._domain

for i,builder in enumerate([a, b, c, d]):
    print(f"domain {'abcd'[i]}: {builder._domain}")


a.node(1, 3.0, 4.0)
b.node(1, 3.0, 4.0)

assert a._domain.getNode(1) != b._domain.getNode(1)

try:
    c.node(1, 3.0, 4.0)
except:
    print("Failed to add node to domain`c`")
try:
    d.node(1, 3.0, 4.0)
except:
    print("Failed to add node to domain`c`")


