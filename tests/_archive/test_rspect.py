import pyg3
import quakeio
import numpy as np


def main(argc:int, argv:list)->int:

    ## create the nodes using constructor: 
    ##        Node(tag, ndof, crd1, crd2)
    ## and then add them to the domain
    rt.eval("""
    model basic 2 2
    node 1   0.0  0.0
    node 2 144.0  0.0
    """)

    # create an elastic material using constriuctor:  
    #    ElasticMaterialModel(tag, E)
    #-- theMaterial = pyg3.ElasticMaterial(1, 3000)

    ## create the truss elements using constructor:
    ##        Truss(tag, dim, nd1, nd2, Material &,A)
    ## and then add them to the domain

    #-- Truss *truss1 = pyg3.Truss(1, 2, 1, 4, *theMaterial, 10.0)

    #-- theDomain.addElement(truss1)
    
    # create the single-point constraint objects using constructor:
    #   SP_Constraint(tag, nodeTag, dofID, value)
    # and then add them to the domain
    
    sp1 = pyg3.domain.SP_Constraint(1, 1, 0, 0.0)
    print(sp1)
    print(theDomain.addSP_Constraint(sp1))
    #sp2 = pyg3.domain.SP_Constraint(2, 1, 1, 0.0)    
    #theDomain.addSP_Constraint(sp2)

    # construct a linear time series object using constructor:
    #   LinearSeries()

    ground_motion = quakeio.read("IELC180.AT2")
    accel = np.asarray(ground_motion) 
    # construct a load pattren using constructor:
    #     LoadPattern(tag)
    # and then set it's TimeSeries and add it to the domain 
    theLoadPattern = pyg3.domain.UniformExcitation(
        tag       = 1,
        dof       = 2,
        accel     = accel,
        time_step = ground_motion.time_step
    )
    print(theLoadPattern)
    theDomain.addLoadPattern(theLoadPattern)
    print(theDomain)


    # create an Analysis object to perform a static analysis of the model
    #  - constructs:
    #    AnalysisModel of type AnalysisModel,
    #      EquiSolnAlgo of type Linear
    #      StaticIntegrator of type LoadControl
    #      ConstraintHandler of type Penalty
    #    DOF_Numberer which uses RCM
    #    LinearSOE of type Band SPD
    # and then the StaticAnalysis object    
    theModel      = pyg3.analysis.AnalysisModel()
    print(theModel)
    theSolnAlgo   = pyg3.analysis.Linear()
    print(theSolnAlgo)
    theIntegrator = pyg3.analysis.LoadControl(1.0, 1, 1.0, 1.0)
    theHandler    = pyg3.analysis.PenaltyConstraintHandler(1.0e8,1.0e8)
    theRCM        = pyg3.graph.RCM()
    print(theRCM)
    theNumberer   = pyg3.analysis.DOF_Numberer(theRCM)
    print(theNumberer)
    theSolver     = pyg3.sys_of_eqn.BandSPDLinLapackSolver()
    print(theSolver)
    theSOE        = pyg3.sys_of_eqn.BandSPDLinSOE(theSolver)

    theAnalysis   = pyg3.analysis.DirectIntegrationAnalysis(
                     theDomain,
                     theHandler,
                     theNumberer,
                     theModel,
                     theSolnAlgo,
                     theSOE,
                     theIntegrator
                   )
    print(theAnalysis)

    # perform the analysis & print out the results for the domain
    numSteps = 1
    theAnalysis.analyze(numSteps)
    print("Analysis complete")
#    return 0

main(0,0)



