# written: fmk
# date: 02/99
#
# purpose: example1 in OpenSeesIntro.tex

# add nodes - command: node nodeId xCrd yCrd
model = {
    "nodes": {
         1: [  0.0,  0.0],
         2: [144.0,  0.0],
         3: [168.0,  0.0],
         4: [ 72.0, 96.0]
    },
    "supports": {
        1:  [1, 1]
        2:  [1, 1]
        3:  [1, 1]
    }
}

# add material - command: uniaxialMaterial <matType> matID <matArgs>
uniaxialMaterial Elastic 1 3000

# add truss elements - command: element truss trussID node1 node2 A matID
element truss 1 1 4 10.0 3000
element truss 2 2 4  5.0 3000
element truss 3 3 4  5.0 3000

#create the ModelBuilder object
model BasicBuilder -ndm 2 -ndf 2

# build the model 
# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?

pattern Plain 1 "Linear" {
    # apply the load - command: load nodeID xForce yForce
    load 4 100 -50
}


# build the components for the analysis object

# create the analysis object 
analysis Static {
    system BandSPD
    constraints Plain
    integrator LoadControl 1.0
    algorithm Linear
    numberer RCM
}

