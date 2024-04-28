import numpy as np
import opensees.openseespy

# Create the model

def arch_model():

    # Define model parameters
    L      = 5000
    Rise   = 500
    Offset = 200

    # Define material parameters
    E = 200
    A = 1e4
    I = 1e8

    # Compute radius
    R  = Rise/2 + (2*L)**2/(8*Rise)
    th = 2*np.arcsin(L/R)

    #
    # Build the model
    #
    model = opensees.openseespy.Model(ndm=2, ndf=3)

    # Create nodes
    ne  = 10
    nen = 2
    nn  = ne*(nen-1)+1
    mid = (nn+1)//2      # midpoint node

    for i, angle in enumerate(np.linspace(-th/2, th/2, nn)):
        tag = i + 1

        # Compute x and add offset if midpoint
        x = R*np.sin(angle)
        if tag == mid:
            x -= Offset

        # Compute y
        y = R*np.cos(angle)

        # create the node
        model.node(tag, x, y)


    # Create elements
    transfTag = 1
    model.geomTransf("Corotational", transfTag)
    for i in range(ne):
        tag   = i+1
        nodes = (i+1, i+2)
        model.element("ElasticBeamColumn", tag, *nodes, A, E, I, transfTag)


    model.fix( 1, 1, 1, 0)
    model.fix(nn, 1, 1, 0)
    
    # Create a load pattern that scales linearly
    model.pattern("Plain", 1, "Linear")

    # Add a nodal load to the pattern
    model.load(mid, 0.0, -1.0, 0.0, pattern=1)

    
    model.system("ProfileSPD")
    # model.system("FullGeneral")
    # model.system("BandGeneral")
    # model.system("Umfpack", det=True)

    model.test("NormUnbalance", 1e-6, 25, 0)
    model.algorithm("Newton")
    model.analysis("Static")


    return model, mid


