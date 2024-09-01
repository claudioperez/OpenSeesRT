import opensees
mat_1 = ...
mat_2 = ...

block_1  = {
    1: (0.0, 0.0),
    2: (1.1, 0.0),
    3: (1.0, 1.0),
    4: (0.0, 1.0),
    5: (0.5,-0.1),
    6: (1.1, 0.5)
}


block_2  = {
    1: (1.1, 0.0),
    2: (2.0, 0.0),
    3: (2.0, 1.0),
    4: (1.0, 1.0),
    5: (1.5,-0.1),
#   7: (2.1, 0.5),
    8: (1.1, 0.5)
}

model = opensees.model(ndm=2, ndf=2)

model.block((2,2), "ShellMITC4", ["1", "PlaneStress", "1"], block_1)

model.block((2,2), "ShellMITC4", ["2", "PlaneStress", "1"], block_2)


#print(opensees.tcl.dumps(model))



