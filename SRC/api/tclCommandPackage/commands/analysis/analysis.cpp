#include <g3_api.h>
#include <runtimeAPI.h>
#include <analysisAPI.h>
#include <OPS_Globals.h>
#include <G3_Logging.h>

#include <Domain.h> // for modal damping

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// analysis model
#include <AnalysisModel.h>
#include <LoadControl.h>
#include <EquiSolnAlgo.h>

// Eigenvalue analysis
#include <EigenSOE.h>
#include <EigenSolver.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <SymArpackSOE.h>
#include <SymArpackSolver.h>
#include <BandArpackSOE.h>
#include <BandArpackSolver.h>
#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <FullGenEigenSolver.h>
// for response spectrum analysis
extern void OPS_DomainModalProperties(G3_Runtime*);
extern void OPS_ResponseSpectrumAnalysis(G3_Runtime*);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char **argv, Domain *domain);
static int numEigen = 0;


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
//  #include <PenaltyHandlerNoHomoSPMultipliers.h>
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

extern EquiSolnAlgo *theAlgorithm ;
extern ConstraintHandler *theHandler ;
extern DOF_Numberer *theNumberer ;
// extern LinearSOE *theSOE ;
extern EigenSOE *theEigenSOE ;
extern TransientIntegrator *theTransientIntegrator;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern VariableTimeStepDirectIntegrationAnalysis
           *theVariableTimeStepTransientAnalysis;
extern ConvergenceTest *theTest;

LinearSOE *G3_getDefaultLinearSoe(G3_Runtime* rt, int flags);

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
  LinearSOE *theSOE = G3_getDefaultLinearSoe(rt, 0);

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
      opswrn << "analysis Static - no Algorithm yet specified, \n"
             << " NewtonRaphson default will be used\n";

      theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
      opswrn
          << "WARNING analysis Static - no ConstraintHandler yet specified, \n"
             << " PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
      opswrn << "analysis Static - no Numberer specified, \n"
             << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }
    if (the_static_integrator == 0) {
      opswrn << "analysis Static - no integrator specified, \n"
             << " StaticIntegrator default will be used\n";
      the_static_integrator = new LoadControl(1, 1, 1, 1);
      G3_setStaticIntegrator(rt, the_static_integrator);
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
      opswrn << "analysis Transient - no Algorithm yet specified, \n"
             << " NewtonRaphson default will be used\n";

      theAlgorithm = new NewtonRaphson(*theTest);
    }
    if (theHandler == 0) {
      opswrn << "analysis Transient dt tFinal - no ConstraintHandler\n"
             << " yet specified, PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }
    if (theNumberer == 0) {
      opswrn
          << "WARNING analysis Transient dt tFinal - no Numberer specified, \n"
             << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }
    if (theTransientIntegrator == 0) {
      opswrn << "analysis Transient dt tFinal - no Integrator specified, \n"
             << " Newmark(.5,.25) default will be used\n";
      theTransientIntegrator = new Newmark(0.5, 0.25);
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
      opswrn << "analysis Transient - no Algorithm yet specified, \n"
             << " NewtonRaphson default will be used\n";
      theAlgorithm = new NewtonRaphson(*theTest);
    }

    if (theHandler == 0) {
      opswrn << "analysis Transient dt tFinal - no ConstraintHandler\n"
             << " yet specified, PlainHandler default will be used\n";
      theHandler = new PlainHandler();
    }

    if (theNumberer == 0) {
      opswrn
          << "analysis Transient dt tFinal - no Numberer specified, \n"
          << " RCM default will be used\n";
      RCM *theRCM = new RCM(false);
      theNumberer = new DOF_Numberer(*theRCM);
    }

    if (theTransientIntegrator == 0) {
      opswrn << "analysis Transient dt tFinal - no Integrator specified, \n"
             << "Newmark(.5,.25) default will be used\n";
      theTransientIntegrator = new Newmark(0.5, 0.25);
    }
/* cmp
 * changed so that G3_setLinearSoe creates default
    if (theSOE == 0) {
      opserr << "WARNING analysis Transient dt tFinal - no LinearSOE "
                "specified, \n";
      opserr << " ProfileSPDLinSOE default will be used\n";
      ProfileSPDLinSolver *theSolver;
      theSolver = new ProfileSPDLinDirectSolver();
    }
*/

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

int
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  /* static */ char *resDataPtr = 0;
  /* static */ int resDataSize = 0;
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);
  AnalysisModel* the_analysis_model = G3_getAnalysisModel(rt);
  DirectIntegrationAnalysis *directIntAnalysis = G3_getTransientAnalysis(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - eigen <type> numModes?\n";
    return TCL_ERROR;
  }

  bool generalizedAlgo =
      true; // 0 - frequency/generalized (default),1 - standard, 2 - buckling
  int typeSolver = EigenSOE_TAGS_ArpackSOE;
  int loc = 1;
  double shift = 0.0;
  bool findSmallest = true;

  // Check type of eigenvalue analysis
  while (loc < (argc - 1)) {
    if ((strcmp(argv[loc], "frequency") == 0) ||
        (strcmp(argv[loc], "-frequency") == 0) ||
        (strcmp(argv[loc], "generalized") == 0) ||
        (strcmp(argv[loc], "-generalized") == 0))
      generalizedAlgo = true;

    else if ((strcmp(argv[loc], "standard") == 0) ||
             (strcmp(argv[loc], "-standard") == 0))
      generalizedAlgo = false;

    else if ((strcmp(argv[loc], "-findLargest") == 0))
      findSmallest = false;

    else if ((strcmp(argv[loc], "genBandArpack") == 0) ||
             (strcmp(argv[loc], "-genBandArpack") == 0) ||
             (strcmp(argv[loc], "genBandArpackEigen") == 0) ||
             (strcmp(argv[loc], "-genBandArpackEigen") == 0))
      typeSolver = EigenSOE_TAGS_ArpackSOE;

    else if ((strcmp(argv[loc], "symmBandLapack") == 0) ||
             (strcmp(argv[loc], "-symmBandLapack") == 0) ||
             (strcmp(argv[loc], "symmBandLapackEigen") == 0) ||
             (strcmp(argv[loc], "-symmBandLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;

    else if ((strcmp(argv[loc], "fullGenLapack") == 0) ||
             (strcmp(argv[loc], "-fullGenLapack") == 0) ||
             (strcmp(argv[loc], "fullGenLapackEigen") == 0) ||
             (strcmp(argv[loc], "-fullGenLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_FullGenEigenSOE;

    else {
      opserr << "eigen - unknown option specified " << argv[loc] << endln;
    }

    loc++;
  }

  // check argv[loc] for number of modes
  if ((Tcl_GetInt(interp, argv[loc], &numEigen) != TCL_OK) || numEigen < 0) {
    opserr << "WARNING eigen numModes?  - illegal numModes\n";
    return TCL_ERROR;
  }

  //
  // create a transient analysis if no analysis exists
  //

  if (the_static_analysis == 0 && directIntAnalysis == 0) {
    if (the_analysis_model == 0)
      the_analysis_model = new AnalysisModel();
    if (theTest == 0)
      theTest = new CTestNormUnbalance(1.0e-6, 25, 0);
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
      theTransientIntegrator = new Newmark(0.5, 0.25);
    }
    LinearSOE *theSOE = G3_getDefaultLinearSoe(rt, 0);
    directIntAnalysis = new DirectIntegrationAnalysis(
        *domain, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
        *theSOE, *theTransientIntegrator, theTest);
  }

  //
  // create a new eigen system and solver
  //
  bool setEigen = false;
  if (theEigenSOE != 0) {
    if (theEigenSOE->getClassTag() != typeSolver) {
      //	delete theEigenSOE;
      theEigenSOE = 0;
      setEigen = true;
    }
  } else {
    if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
      SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
      theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *the_analysis_model);

    } else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {

      FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
      theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *the_analysis_model);

    } else {

      theEigenSOE = new ArpackSOE(shift);
    }

    //
    // set the eigen soe in the system
    //

    if (the_static_analysis != 0) {
      the_static_analysis->setEigenSOE(*theEigenSOE);
    } else if (directIntAnalysis != 0) {
      directIntAnalysis->setEigenSOE(*theEigenSOE);
    }

#ifdef _PARALLEL_PROCESSING

    if (OPS_PARTITIONED == false && OPS_NUM_SUBDOMAINS > 1) {
      if (partitionModel(0) < 0) {
        opserr
            << "WARNING before analysis; partition failed - too few elements\n";
        OpenSeesExit(clientData, interp, argc, argv);
        return TCL_ERROR;
      }
    }

    if (the_static_analysis != 0 || directIntAnalysis != 0) {
      SubdomainIter &theSubdomains = domain->getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
        theSub->setAnalysisEigenSOE(*theEigenSOE);
      }
    }
#endif

  } // theEigenSOE != 0

  int requiredDataSize = 40 * numEigen;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != 0) {
      delete[] resDataPtr;
    }
    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }

  for (int i = 0; i < requiredDataSize; i++)
    resDataPtr[i] = '\n';

  int result = 0;

  if (the_static_analysis != 0) {
    result = the_static_analysis->eigen(numEigen, generalizedAlgo, findSmallest);
  } else if (directIntAnalysis != 0) {
    result =
        directIntAnalysis->eigen(numEigen, generalizedAlgo, findSmallest);
  }

  if (result == 0) {
    //      char *eigenvalueS = new char[15 * numEigen];
    const Vector &eigenvalues = domain->getEigenvalues();
    int cnt = 0;
    for (int i = 0; i < numEigen; i++) {
      cnt += sprintf(&resDataPtr[cnt], "%35.20f  ", eigenvalues[i]);
    }

    Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
  }

  return TCL_OK;
}

int
modalProperties(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, the_domain);
  OPS_DomainModalProperties(rt);
  return TCL_OK;
}

int
responseSpectrum(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, the_domain);
  OPS_ResponseSpectrumAnalysis(rt);
  return TCL_OK;
}

int
modalDamping(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  if (argc < 2) {
    opserr
        << "WARNING modalDamping ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  if (numEigen == 0 || theEigenSOE == 0) {
    opserr << "WARNING - modalDmping - eigen command needs to be called first "
              "- NO MODAL DAMPING APPLIED\n ";
  }

  int numModes = argc - 1;
  double factor;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    opserr << "WARNING modalDmping - same #damping factors as modes must be "
              "specified\n";
    opserr
        << "                    - same damping ratio will be applied to all\n";
  }

  //
  // read in values and set factors
  //

  if (numModes == numEigen) {

    for (int i = 0; i < numEigen; i++) {
      if (Tcl_GetDouble(interp, argv[1 + i], &factor) != TCL_OK) {
        opserr << "WARNING modalDamping - could not read factor for model "
               << i + 1 << endln;
        return TCL_ERROR;
      }
      modalDampingValues[i] = factor;
    }

  } else {

    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << "WARNING modalDamping - could not read factor for all modes \n";
      return TCL_ERROR;
    }

    for (int i = 0; i < numEigen; i++)
      modalDampingValues[i] = factor;
  }

  // set factors in domain
  Domain *theDomain = G3_getDomain(G3_getRuntime(interp));
  theDomain->setModalDampingFactors(&modalDampingValues, true);

  // opserr << "modalDamping Factors: " << modalDampingValues;

  return TCL_OK;
}

int
modalDampingQ(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  if (argc < 2) {
    opserr
        << "WARNING modalDamping ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  if (numEigen == 0 || theEigenSOE == 0) {
    opserr << "WARINING - modalDmping - eigen command needs to be called first "
              "- NO MODAL DAMPING APPLIED\n ";
  }

  int numModes = argc - 1;
  double factor = 0;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    opserr << "WARNING modalDmping - same #damping factors as modes must be "
              "specified\n";
    opserr << "                    - same damping ratio will be applied to all";
  }

  //
  // read in values and set factors
  //

  if (numModes == numEigen) {

    // read in all factors one at a time
    for (int i = 0; i < numEigen; i++) {
      if (Tcl_GetDouble(interp, argv[1 + i], &factor) != TCL_OK) {
        opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
                  "read betaK? \n";
        return TCL_ERROR;
      }
      modalDampingValues[i] = factor;
    }

  } else {

    //  read in one & set all factors to that value
    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
                "read betaK? \n";
      return TCL_ERROR;
    }

    for (int i = 0; i < numEigen; i++)
      modalDampingValues[i] = factor;
  }

  // set factors in domain
  Domain *theDomain = G3_getDomain(G3_getRuntime(interp));
  theDomain->setModalDampingFactors(&modalDampingValues, false);
  return TCL_OK;
}

