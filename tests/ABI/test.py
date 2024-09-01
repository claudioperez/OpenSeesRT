import pyg3
exec(open("twin").read())

a = pyg3.AbstractAnalysis(structure.model)

a.build(k_tran=3.0, k_sect=390.0, k_vert=390, k_long=30.0, k_plan=30.0, k_elev=390.0)

a.recorder("Node", "-nodeRange", 0, len(structure.model.nodes), "eigen 1")

a.eigen(1)
#a.nodeEigenvector(0, 1)

domain = pyg3._pyg3.get_domain(a._interp.tk.interpaddr())

n = domain.getNode(1)
v = n.getDisp()

pyg3._pyg3.copy_vector(v)

print(n.getEigenvectors())

del n
print("node deleted")
del domain

del a

