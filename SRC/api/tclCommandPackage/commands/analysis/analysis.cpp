
#include <g3_api.h>
#include <runtimeAPI.h>
#include <analysisAPI.h>
#include <OPS_Globals.h>

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// analysis model
#include <AnalysisModel.h>
#include <LoadControl.h>
#include <EquiSolnAlgo.h>

// integrators
#include <Newmark.h>
#include <StagedNewmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Newmark1.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>

#include <LoadControl.h>
#include <StagedLoadControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
#include <HSConstraint.h>
#include <MinUnbalDispNorm.h>
#include <DisplacementControl.h>
#include <EQPath.h>


// convergence tests
#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>
#include <CTestRelativeNormUnbalance.h>
#include <CTestRelativeNormDispIncr.h>
#include <CTestRelativeEnergyIncr.h>
#include <CTestRelativeTotalNormDispIncr.h>
#include <CTestFixedNumIter.h>
#include <NormDispAndUnbalance.h>
#include <NormDispOrUnbalance.h>
#ifdef OPS_USE_PFEM
#  include <CTestPFEM.h>
#endif

// soln algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <NewtonLineSearch.h>
#include <ModifiedNewton.h>
#include <Broyden.h>
#include <BFGS.h>
#include <KrylovNewton.h>
#include <PeriodicNewton.h>
#include <AcceleratedNewton.h>
#include <ExpressNewton.h>

// constraint handlers
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
//#include <PenaltyHandlerNoHomoSPMultipliers.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <PlainNumberer.h>
#include <DOF_Numberer.h>
// graph
#include <RCM.h>
#include <AMDNumberer.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

#include <ConjugateGradientSolver.h>

#ifdef _ITPACK
// #include <ItpackLinSOE.h>
// #include <ItpackLinSolver.h>
#endif

#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>

#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DiagonalSOE.h>
#include <DiagonalDirectSolver.h>

#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>



extern EquiSolnAlgo *theAlgorithm ;
extern ConstraintHandler *theHandler ;
extern DOF_Numberer *theNumberer ;
extern LinearSOE *theSOE ;
extern EigenSOE *theEigenSOE ;
extern TransientIntegrator *theTransientIntegrator;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern VariableTimeStepDirectIntegrationAnalysis
           *theVariableTimeStepTransientAnalysis;
extern ConvergenceTest *theTest;

//
// command invoked to allow the Analysis object to be built
int
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);

  StaticIntegrator *the_static_integrator = G3_getStaticIntegrator(rt);
  AnalysisModel* the_analysis_model = nullptr;

  // make sure at least one other argument to contain type of system

  if (argc < 2) {
    opserr << "WARNING need to specify an analysis type (Static, Transient)\n";
    return TCL_ERROR;
  }

  // do nothing if request is for the same analysis type!

  if ((strcmp(argv[1], "Static") == 0) && (the_static_analysis != 0))
    return TCL_OK;

  if (((strcmp(argv[1], "VariableTimeStepTransient") == 0) ||
       (strcmp(argv[1], "TransientWithVariableTimeStep") == 0) ||
       (strcmp(argv[1], "VariableTransient") == 0)) &&
      (theVariableTimeStepTransientAnalysis != 0))
    return TCL_OK;

  if ((strcmp(argv[1], "Transient") == 0) && (theTransientAnalysis != 0))
    return TCL_OK;

  // analysis changing .. delete the old analysis

  if (the_static_analysis != 0) {
    G3_delStaticAnalysis(rt);
    delete the_static_analysis;
    the_static_analysis = 0;
    opserr << "WARNING: analysis .. StaticAnalysis already exists => "
              "wipeAnalysis not invoked, problems may arise\n";
  }

  if (theTransientAnalysis != 0) {
    delete theTransientAnalysis;
    theTransientAnalysis = 0;
    theVariableTimeStepTransientAnalysis = 0;
    opserr << "WARNING: analysis .. TransientAnalysis already exists => "
              "wipeAnalysis not invoked, problems may arise\n";
  }

  // check argv[1] for type of SOE and create it
  if (strcmp(argv[1], "Static") == 0) {
    the_analysis_model = G3_getAnalysisModel(rt);
    // make sure all the components have been built,
    // otherwise print a warning and use some defaults
    if (the_analysis_model == 0){
      the_analysis_model = new AnalysisModel();
      G3_setAnalysisModel(rt, the_analysis_model);
    }

    if (theTest == 0)
      theTest = new CTestNormUnbalance(1.0e-6, 25, 0);

    if (theAlgorithm == 0) {
      opserr << "WARNING analysis Static - no Algorithm yet specified, \n";
      opserr << " NewtonRaphson default will be used\n";

      theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
      opserr
          << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
      opserr << " PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
      // opserr << "WARNING analysis Static - no Numberer specified, \n";
      // opserr << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }
    if (the_static_integrator == 0) {
      opserr << "WARNING analysis Static - no integrator specified, \n";
      opserr << " StaticIntegrator default will be used\n";
      the_static_integrator = new LoadControl(1, 1, 1, 1);
      G3_setStaticIntegrator(rt, the_static_integrator);
    }
    if (theSOE == 0) {
      opserr << "WARNING analysis Static - no LinearSOE specified, \n";
      opserr << " ProfileSPDLinSOE default will be used\n";
      ProfileSPDLinSolver *theSolver;
      theSolver = new ProfileSPDLinDirectSolver();
#ifdef _PARALLEL_PROCESSING
      theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
      theSOE = new ProfileSPDLinSOE(*theSolver);
#endif
    }

    the_static_analysis = new StaticAnalysis(
        *domain, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
        *theSOE, *the_static_integrator, theTest);

    G3_setStaticAnalysis(rt, the_static_analysis);


// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
    if (theSensitivityAlgorithm != 0 &&
        theSensitivityAlgorithm->shouldComputeAtEachStep()) {
      // the_static_analysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
    }
#endif
    // AddingSensitivity:END /////////////////////////////////
#ifdef OPS_USE_PFEM
  } else if (strcmp(argv[1], "PFEM") == 0) {

    if (argc < 5) {
      opserr << "WARNING: wrong no of args -- analysis PFEM dtmax dtmin "
                "gravity <ratio>\n";
      return TCL_ERROR;
    }
    double dtmax, dtmin, gravity, ratio = 0.5;
    if (Tcl_GetDouble(interp, argv[2], &dtmax) != TCL_OK) {
      opserr << "WARNING: invalid dtmax " << argv[2] << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &dtmin) != TCL_OK) {
      opserr << "WARNING: invalid dtmin " << argv[3] << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &gravity) != TCL_OK) {
      opserr << "WARNING: invalid gravity " << argv[4] << "\n";
      return TCL_ERROR;
    }
    if (argc > 5) {
      if (Tcl_GetDouble(interp, argv[5], &ratio) != TCL_OK) {
        opserr << "WARNING: invalid ratio " << argv[5] << "\n";
        return TCL_ERROR;
      }
    }

    if (the_analysis_model == 0) {
      the_analysis_model = new AnalysisModel();
      G3_setAnalysisModel(rt, the_analysis_model);
    }
    if (theTest == 0) {
      // theTest = new CTestNormUnbalance(1e-2,10000,1,2,3);
      theTest =
          new CTestPFEM(1e-2, 1e-2, 1e-2, 1e-2, 1e-4, 1e-3, 10000, 100, 1, 2);
    }
    if (theAlgorithm == 0) {
      theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
      theHandler = new TransformationConstraintHandler();
    }
    if (theNumberer == 0) {
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }
    if (theTransientIntegrator == 0) {
      theTransientIntegrator = new PFEMIntegrator();
    }
    if (theSOE == 0) {
      PFEMSolver *theSolver = new PFEMSolver();
      theSOE = new PFEMLinSOE(*theSolver);
    }
    thePFEMAnalysis = new PFEMAnalysis(theDomain, *theHandler, *theNumberer,
                                       *the_analysis_model, *theAlgorithm,
                                       *theSOE, *theTransientIntegrator,
                                       theTest, dtmax, dtmin, gravity, ratio);

    theTransientAnalysis = thePFEMAnalysis;
#endif
  } else if (strcmp(argv[1], "Transient") == 0) {
    // make sure all the components have been built,
    // otherwise print a warning and use some defaults
    if (the_analysis_model == 0){
      the_analysis_model = new AnalysisModel();
      G3_setAnalysisModel(rt, the_analysis_model);
    }

    if (theTest == 0)
      theTest = new CTestNormUnbalance(1.0e-6, 25, 0);

    if (theAlgorithm == 0) {
      opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
      opserr << " NewtonRaphson default will be used\n";

      theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
      opserr << " yet specified, PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
      opserr
          << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
      opserr << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }
    if (theTransientIntegrator == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no Integrator "
                "specified, \n";
      opserr << " Newmark(.5,.25) default will be used\n";
      theTransientIntegrator = new Newmark(0.5, 0.25);
    }
    if (theSOE == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no LinearSOE "
                "specified, \n";
      opserr << " ProfileSPDLinSOE default will be used\n";
      ProfileSPDLinSolver *theSolver;
      theSolver = new ProfileSPDLinDirectSolver();
    }

    int count = 2;
    int numSubLevels = 0;
    int numSubSteps = 10;
    while (count < argc) {
      if (strcmp(argv[count], "-numSubLevels") == 0) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &numSubLevels) != TCL_OK)
            return TCL_ERROR;
      } else if ((strcmp(argv[count], "-numSubSteps") == 0)) {
        count++;
        if (count < argc)
          if (Tcl_GetInt(interp, argv[count], &numSubSteps) != TCL_OK)
            return TCL_ERROR;
      }
      count++;
    }

    theTransientAnalysis = new DirectIntegrationAnalysis(
        *domain, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
        *theSOE, *theTransientIntegrator, theTest, numSubLevels, numSubSteps);
    ;

// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
    if (theSensitivityAlgorithm != 0 &&
        theSensitivityAlgorithm->shouldComputeAtEachStep()) {

      /* This if-statement cannot possibly stay in the code -- MHS
      if(theSensitivityAlgorithm->newAlgorithm()){
        opserr << "WARNING original sensitivity algorothm needs to be specified
      \n"; opserr << "for static analysis \n"; return TCL_ERROR;
      }
      */

      // theTransientAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
    }
#endif
    // AddingSensitivity:END /////////////////////////////////

  } else if ((strcmp(argv[1], "VariableTimeStepTransient") == 0) ||
             (strcmp(argv[1], "TransientWithVariableTimeStep") == 0) ||
             (strcmp(argv[1], "VariableTransient") == 0)) {
    // make sure all the components have been built,
    // otherwise print a warning and use some defaults
    if (the_analysis_model == 0)
      the_analysis_model = new AnalysisModel();

    if (theTest == 0)
      theTest = new CTestNormUnbalance(1.0e-6, 25, 0);

    if (theAlgorithm == 0) {
      opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
      opserr << " NewtonRaphson default will be used\n";
      theAlgorithm = new NewtonRaphson(*theTest);
    }

    if (theHandler == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
      opserr << " yet specified, PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }

    if (theNumberer == 0) {
      opserr
          << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
      opserr << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }

    if (theTransientIntegrator == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no Integrator "
                "specified, \n";
      opserr << " Newmark(.5,.25) default will be used\n";
      theTransientIntegrator = new Newmark(0.5, 0.25);
    }

    if (theSOE == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no LinearSOE "
                "specified, \n";
      opserr << " ProfileSPDLinSOE default will be used\n";
      ProfileSPDLinSolver *theSolver;
      theSolver = new ProfileSPDLinDirectSolver();
    }

    theVariableTimeStepTransientAnalysis =
        new VariableTimeStepDirectIntegrationAnalysis(
            *domain, *theHandler, *theNumberer, *the_analysis_model,
            *theAlgorithm, *theSOE, *theTransientIntegrator, theTest);

    // set the pointer for variabble time step analysis
    theTransientAnalysis = theVariableTimeStepTransientAnalysis;

  } else {
    opserr << "WARNING No Analysis type exists (Static Transient only) \n";
    return TCL_ERROR;
  }

  if (theEigenSOE != 0) {
    if (the_static_analysis != 0 ) {
      the_static_analysis->setEigenSOE(*theEigenSOE);
    } else if (theTransientAnalysis != 0) {
      theTransientAnalysis->setEigenSOE(*theEigenSOE);
    }
  }
  return TCL_OK;
}
