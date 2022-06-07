#include <g3_api.h>
#include <runtimeAPI.h>
#include <analysisAPI.h>
#include <G3Parse.h>

#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

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

// line search
#include <BisectionLineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <SecantLineSearch.h>

// accelerators
#include <RaphsonAccelerator.h>
#include <PeriodicAccelerator.h>
#include <KrylovAccelerator.h>
#include <SecantAccelerator1.h>
#include <SecantAccelerator2.h>
#include <SecantAccelerator3.h>
//#include <MillerAccelerator.h>


extern EquiSolnAlgo *theAlgorithm;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern StaticAnalysis *theStaticAnalysis;
extern ConvergenceTest *theTest;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char **argv, Domain *domain);


EquiSolnAlgo*
G3Parse_newEquiSolnAlgo(G3_Builder* rt, int argc, G3_Char** argv);
EquiSolnAlgo*
G3Parse_newSecantNewtonAlgorithm(G3_Builder*, int, G3_Char **);
EquiSolnAlgo*
G3Parse_newLinearAlgorithm(G3_Builder*, int, G3_Char **);

//
// command invoked to allow the SolnAlgorithm object to be built
//
int
specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{

  // make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << "WARNING need to specify an Algorithm type \n";
    return TCL_ERROR;
  }
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);
  EquiSolnAlgo *theNewAlgo = 0;
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, domain);

  theNewAlgo = G3Parse_newEquiSolnAlgo(rt, argc, argv);
  if (theNewAlgo == nullptr) {
    return TCL_ERROR;

  } else {
    if (theTest != 0)
      theNewAlgo->setConvergenceTest(theTest);
    theAlgorithm = theNewAlgo;

    // if the analysis exists - we want to change the SOE
    if (the_static_analysis != 0)
      the_static_analysis->setAlgorithm(*theAlgorithm);
    else if (theTransientAnalysis != 0)
      theTransientAnalysis->setAlgorithm(*theAlgorithm);

#ifdef _PARALLEL_PROCESSING
    if (the_static_analysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
        theSub->setAnalysisAlgorithm(*theAlgorithm);
      }
    }
#endif
  }

  return TCL_OK;
}

EquiSolnAlgo*
G3Parse_newEquiSolnAlgo(G3_Builder* rt, int argc, G3_Char** argv)
{
  EquiSolnAlgo *theNewAlgo = nullptr;

  // check argv[1] for type of Algorithm and create the object
  if (strcmp(argv[1], "Linear") == 0) {
    theNewAlgo = G3Parse_newLinearAlgorithm(rt, argc, argv);
  }

  else if (strcmp(argv[1], "Newton") == 0) {
    void *theNewtonAlgo = OPS_NewtonRaphsonAlgorithm(rt);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if ((strcmp(argv[1], "NewtonHallM") == 0) ||
           (strcmp(argv[1], "NewtonHall") == 0)) {
    void *theNewtonAlgo = OPS_NewtonHallM(rt);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if (strcmp(argv[1], "ModifiedNewton") == 0) {
    void *theNewtonAlgo = OPS_ModifiedNewton(rt);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else if (strcmp(argv[1], "SecantNewton") == 0) {
    theNewAlgo = G3Parse_newSecantNewtonAlgorithm(rt, argc, argv);
  }

/*
  else if (strcmp(argv[1], "KrylovNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          iterateTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          iterateTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          iterateTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          incrementTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          incrementTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          incrementTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
        i++;
        maxDim = atoi(argv[i]);
      }
    }
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    Accelerator *theAccel;
    theAccel = new KrylovAccelerator(maxDim, iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1], "RaphsonNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          iterateTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          iterateTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          iterateTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          incrementTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          incrementTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          incrementTangent = NO_TANGENT;
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    Accelerator *theAccel;
    theAccel = new RaphsonAccelerator(iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1], "MillerNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;

    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          iterateTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          iterateTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          iterateTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          incrementTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          incrementTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          incrementTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
        i++;
        maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    Accelerator *theAccel = 0;
    // theAccel = new MillerAccelerator(maxDim, 0.01, iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1], "PeriodicNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          iterateTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          iterateTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          iterateTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
        i++;
        if (strcmp(argv[i], "current") == 0)
          incrementTangent = CURRENT_TANGENT;
        if (strcmp(argv[i], "initial") == 0)
          incrementTangent = INITIAL_TANGENT;
        if (strcmp(argv[i], "noTangent") == 0)
          incrementTangent = NO_TANGENT;
      } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
        i++;
        maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    Accelerator *theAccel;
    theAccel = new PeriodicAccelerator(maxDim, iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1], "Broyden") == 0) {
    int formTangent = CURRENT_TANGENT;
    int count = -1;

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-secant") == 0) {
        formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i], "-initial") == 0) {
        formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
        count = atoi(argv[i]);
      }
    }

    if (count == -1)
      theNewAlgo = new Broyden(*theTest, formTangent);
    else
      theNewAlgo = new Broyden(*theTest, formTangent, count);
  }

  else if (strcmp(argv[1], "BFGS") == 0) {
    int formTangent = CURRENT_TANGENT;
    int count = -1;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-secant") == 0) {
        formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i], "-initial") == 0) {
        formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++], "-count") == 0 && i < argc) {
        count = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    if (count == -1)
      theNewAlgo = new BFGS(*theTest, formTangent);
    else
      theNewAlgo = new BFGS(*theTest, formTangent, count);
  }

  else if (strcmp(argv[1], "NewtonLineSearch") == 0) {
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return nullptr;
    }

    int count = 2;

    // set some default variable
    double tol = 0.8;
    int maxIter = 10;
    double maxEta = 10.0;
    double minEta = 0.1;
    int pFlag = 1;
    int typeSearch = 0;

    while (count < argc) {
      if (strcmp(argv[count], "-tol") == 0) {
        count++;
        if (Tcl_GetDouble(interp, argv[count], &tol) != TCL_OK)
          return nullptr;
        count++;
      } else if (strcmp(argv[count], "-maxIter") == 0) {
        count++;
        if (G3Parse_getInt(rt, argv[count], &maxIter) != TCL_OK)
          return nullptr;
        count++;
      } else if (strcmp(argv[count], "-pFlag") == 0) {
        count++;
        if (G3Parse_getInt(rt, argv[count], &pFlag) != TCL_OK)
          return nullptr;
        count++;
      } else if (strcmp(argv[count], "-minEta") == 0) {
        count++;
        if (Tcl_GetDouble(interp, argv[count], &minEta) != TCL_OK)
          return nullptr;
        count++;
      } else if (strcmp(argv[count], "-maxEta") == 0) {
        count++;
        if (Tcl_GetDouble(interp, argv[count], &maxEta) != TCL_OK)
          return nullptr;
        count++;
      } else if (strcmp(argv[count], "-type") == 0) {
        count++;
        if (strcmp(argv[count], "Bisection") == 0)
          typeSearch = 1;
        else if (strcmp(argv[count], "Secant") == 0)
          typeSearch = 2;
        else if (strcmp(argv[count], "RegulaFalsi") == 0)
          typeSearch = 3;
        else if (strcmp(argv[count], "LinearInterpolated") == 0)
          typeSearch = 3;
        else if (strcmp(argv[count], "InitialInterpolated") == 0)
          typeSearch = 0;
        count++;
      } else
        count++;
    }

    LineSearch *theLineSearch = 0;
    if (typeSearch == 0)
      theLineSearch = new InitialInterpolatedLineSearch(tol, maxIter, minEta,
                                                        maxEta, pFlag);

    else if (typeSearch == 1)
      theLineSearch =
          new BisectionLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 2)
      theLineSearch = new SecantLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 3)
      theLineSearch =
          new RegulaFalsiLineSearch(tol, maxIter, minEta, maxEta, pFlag);

    theNewAlgo = new NewtonLineSearch(*theTest, theLineSearch);
  }
*/

  else if (strcmp(argv[1], "ExpressNewton") == 0) {
    void *theNewtonAlgo = OPS_ExpressNewton(rt);
    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
  }

  else {
    opserr << "WARNING No EquiSolnAlgo type " << argv[1] << " exists\n";
    return nullptr;
  }

  return theNewAlgo;
}

int
sensitivityIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                      TCL_Char **argv)
{
  // Does nothing, but keeping command for backward compatibility
  return TCL_OK;
}

EquiSolnAlgo*
G3Parse_newLinearAlgorithm(G3_Builder* builder, int argc, G3_Char **argv)
{
    int formTangent = CURRENT_TANGENT;
    int factorOnce = 0;
    int count = 2;
    while (count < argc) {
      if ((strcmp(argv[count], "-secant") == 0) ||
          (strcmp(argv[count], "-Secant") == 0)) {
        formTangent = CURRENT_SECANT;
      } else if ((strcmp(argv[count], "-initial") == 0) ||
                 (strcmp(argv[count], "-Initial") == 0)) {
        formTangent = INITIAL_TANGENT;
      } else if ((strcmp(argv[count], "-factorOnce") == 0) ||
                 (strcmp(argv[count], "-FactorOnce") == 0)) {
        factorOnce = 1;
      }
      count++;
    }
    return new Linear(formTangent, factorOnce);
}

EquiSolnAlgo*
G3Parse_newSecantNewtonAlgorithm(G3_Builder* builder, int argc, G3_Char **argv)
{

  ConvergenceTest* theTest = builder->m_global_strategy.m_convergence_test;
  int incrementTangent = CURRENT_TANGENT;
  int iterateTangent = CURRENT_TANGENT;
  int maxDim = 3;
  for (int i = 2; i < argc; i++) {
    if (strcmp(argv[i], "-iterate") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        iterateTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        iterateTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        iterateTangent = NO_TANGENT;
    } else if (strcmp(argv[i], "-increment") == 0 && i + 1 < argc) {
      i++;
      if (strcmp(argv[i], "current") == 0)
        incrementTangent = CURRENT_TANGENT;
      if (strcmp(argv[i], "initial") == 0)
        incrementTangent = INITIAL_TANGENT;
      if (strcmp(argv[i], "noTangent") == 0)
        incrementTangent = NO_TANGENT;
    } else if (strcmp(argv[i], "-maxDim") == 0 && i + 1 < argc) {
      i++;
      maxDim = atoi(argv[i]);
    }
  }

  if (theTest == 0) {
    opserr << "ERROR: No ConvergenceTest yet specified\n";
    return nullptr;
  }

  Accelerator *theAccel;
  theAccel = new SecantAccelerator2(maxDim, iterateTangent);
  return new AcceleratedNewton(*theTest, theAccel, incrementTangent);
}