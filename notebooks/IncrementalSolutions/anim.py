from snap import snap
from sees import render

nstep = 8

model,_ = snap()
mesh = model.asdict()

model.integrator("LoadControl", 400.0)
model.analysis("Static")


for i in range(nstep):
    model.analyze(1)

    resp = {
        i: model.nodeDisp(i) for i in model.getNodeTags()
    }

    render(mesh, resp, noshow=False)


