#include <g3_api.h>
#include <runtimeAPI.h>
#include <analysisAPI.h>
#include <tcl_packages.h>

#include <Domain.h>
#include <Node.h>

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// integrators
#include <LoadControl.h>
#include <StagedLoadControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
#include <HSConstraint.h>
#include <MinUnbalDispNorm.h>
#include <DisplacementControl.h>
#include <EQPath.h>

#include <Newmark.h>
#include <StagedNewmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Newmark1.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>

extern TransientIntegrator *theTransientIntegrator;
extern StaticAnalysis *theStaticAnalysis;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern VariableTimeStepDirectIntegrationAnalysis
           *theVariableTimeStepTransientAnalysis;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char **argv, Domain *domain);

typedef struct externalClassFunction {
  char *funcName;
  void *(*funcPtr)();
  struct externalClassFunction *next;
} ExternalClassFunction;

static ExternalClassFunction *theExternalStaticIntegratorCommands = NULL;
static ExternalClassFunction *theExternalTransientIntegratorCommands = NULL;
static ExternalClassFunction *theExternalAlgorithmCommands = NULL;


//
// command invoked to allow the Integrator object to be built
//
int
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv)
{
  bool assign_to_static_analysis = false;
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);
  StaticIntegrator* the_static_integrator = G3_getStaticIntegrator(rt);
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, domain);

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify an Integrator type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1], "LoadControl") == 0) {
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator LoadControl dlam <Jd "
                "dlamMin dlamMax>\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)
      return TCL_ERROR;
    if (argc > 5) {
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)
        return TCL_ERROR;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    the_static_integrator = new LoadControl(dLambda, numIter, minIncr, maxIncr);

    assign_to_static_analysis = true;
    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  } else if (strcmp(argv[1], "StagedLoadControl") == 0) {
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator StagedLoadControl dlam "
                "<Jd dlamMin dlamMax>\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)
      return TCL_ERROR;
    if (argc > 5) {
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)
        return TCL_ERROR;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    the_static_integrator =
        new StagedLoadControl(dLambda, numIter, minIncr, maxIncr);
    assign_to_static_analysis = true;
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

  else if (strcmp(argv[1], "ArcLength") == 0) {
    double arcLength;
    double alpha;
    if (argc != 4) {
      opserr << "WARNING integrator ArcLength arcLength alpha \n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)
      return TCL_ERROR;
    the_static_integrator = new ArcLength(arcLength, alpha);

    assign_to_static_analysis = true;
    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

  else if (strcmp(argv[1], "ArcLength1") == 0) {
    double arcLength;
    double alpha;
    if (argc != 4) {
      opserr << "WARNING integrator ArcLength1 arcLength alpha \n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)
      return TCL_ERROR;
    the_static_integrator = new ArcLength1(arcLength, alpha);

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }
  /* ************ added for HSConstraint *******************/

  else if (strcmp(argv[1], "HSConstraint") == 0) {
    double arcLength, psi_u, psi_f, u_ref;

    if (argc < 3) {
      opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> "
                "<u_ref> \n";
      return TCL_ERROR;
    }
    if (argc >= 3 && Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
      return TCL_ERROR;
    if (argc >= 4 && Tcl_GetDouble(interp, argv[3], &psi_u) != TCL_OK)
      return TCL_ERROR;
    if (argc >= 5 && Tcl_GetDouble(interp, argv[4], &psi_f) != TCL_OK)
      return TCL_ERROR;
    if (argc == 6 && Tcl_GetDouble(interp, argv[5], &u_ref) != TCL_OK)
      return TCL_ERROR;

    switch (argc) {
    case 3:
      the_static_integrator = new HSConstraint(arcLength);
    case 4:
      the_static_integrator = new HSConstraint(arcLength, psi_u);
    case 5:
      the_static_integrator = new HSConstraint(arcLength, psi_u, psi_f);
    case 6:
      the_static_integrator = new HSConstraint(arcLength, psi_u, psi_f, u_ref);
    }
    assign_to_static_analysis = true;
    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }
  /*********************************************************************************/

  else if (strcmp(argv[1], "MinUnbalDispNorm") == 0) {
    double lambda11, minlambda, maxlambda;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j "
                "maxLambda1j>\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &lambda11) != TCL_OK)
      return TCL_ERROR;
    if (argc > 5) {
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[4], &minlambda) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[5], &maxlambda) != TCL_OK)
        return TCL_ERROR;
    } else {
      minlambda = lambda11;
      maxlambda = lambda11;
      numIter = 1;
      argc += 3;
    }

    int signFirstStepMethod = SIGN_LAST_STEP;
    if (argc == 7)
      if ((strcmp(argv[argc - 1], "-determinant") == 0) ||
          (strcmp(argv[argc - 1], "-det") == 0))
        signFirstStepMethod = CHANGE_DETERMINANT;

    the_static_integrator = new MinUnbalDispNorm(lambda11, numIter, minlambda,
                                               maxlambda, signFirstStepMethod);

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

  else if (strcmp(argv[1], "EQPath") == 0) {
    double arcLength;
    int type;
    int numIter;
    if (argc != 4) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << "REFS : \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
      return TCL_ERROR;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK) {
      opserr << "WARNING integrator $arc_length $type \n";
      opserr << "$type = 1 Minimum Residual Displacement \n";
      opserr << "$type = 2 Normal Plain \n";
      opserr << "$type = 3 Update Normal Plain \n";
      opserr << "$type = 4 Cylindrical Arc-Length \n";

      return TCL_ERROR;
    }

    the_static_integrator = new EQPath(arcLength, type);

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

  else if (strcmp(argv[1], "DisplacementControl") == 0) {
    int node, dof, numIter;
    double increment, minIncr, maxIncr;

    if (argc < 5) {
      opserr << "WARNING integrator DisplacementControl node dof dU \n";
      opserr << "<Jd minIncrement maxIncrement>\n";
      return TCL_ERROR;
    }
    int tangFlag = 0;

    if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)
      return TCL_ERROR;

    if (argc == 6 || argc == 9)
      if (argc == 6) {
        if (strcmp(argv[5], "-initial") == 0)
          tangFlag = 1;
      } else if (strcmp(argv[8], "-initial") == 0)
        tangFlag = 1;

    if (argc > 6) {
      if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK)
        return TCL_ERROR;
    } else {
      minIncr = increment;
      maxIncr = increment;
      numIter = 1;
    }

#ifdef _PARALLEL_PROCESSING

    the_static_integrator = new DistributedDisplacementControl(
        node, dof - 1, increment, numIter, minIncr, maxIncr);
#else
    Node *theNode = domain->getNode(node);
    if (theNode == 0) {
      opserr << "WARNING integrator DisplacementControl node dof dU : Node "
                "does not exist\n";
      return TCL_ERROR;
    }

    int numDOF = theNode->getNumberDOF();
    if (dof <= 0 || dof > numDOF) {
      opserr << "WARNING integrator DisplacementControl node dof dU : invalid "
                "dof given\n";
      return TCL_ERROR;
    }

    the_static_integrator =
        new DisplacementControl(node, dof-1, increment, domain, numIter,
                                minIncr, maxIncr, tangFlag);

    G3_setStaticIntegrator(rt, the_static_integrator);
#endif

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

#ifdef _PARALLEL_INTERPRETERS

  else if ((strcmp(argv[1], "ParallelDisplacementControl") == 0) ||
           (strcmp(argv[1], "ParallelDisplacementControl") == 0)) {
    int node;
    int dof;
    double increment, minIncr, maxIncr;
    int numIter;
    if (argc < 5) {
      opserr << "WARNING integrator DisplacementControl node dof dU \n";
      opserr << "<Jd minIncrement maxIncrement>\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)
      return TCL_ERROR;
    if (argc > 7) {
      if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK)
        return TCL_ERROR;
    } else {
      minIncr = increment;
      maxIncr = increment;
      numIter = 1;
    }

    DistributedDisplacementControl *theDDC = new DistributedDisplacementControl(
        node, dof - 1, increment, numIter, minIncr, maxIncr);

    theDDC->setProcessID(OPS_rank);
    theDDC->setChannels(numChannels, theChannels);
    the_static_integrator = theDDC;

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }
#endif

  else if ((strcmp(argv[1], "TRBDF2") == 0) ||
           (strcmp(argv[1], "Bathe") == 0)) {
    theTransientIntegrator = new TRBDF2();
  }

  else if ((strcmp(argv[1], "TRBDF3") == 0) ||
           (strcmp(argv[1], "Bathe3") == 0)) {
    theTransientIntegrator = new TRBDF3();
  }

  else if (strcmp(argv[1], "Houbolt") == 0) {
    theTransientIntegrator = new Houbolt();
  }

  /*else if (strcmp(argv[1],"ParkLMS3") == 0) {
      theTransientIntegrator = new ParkLMS3();
  }*/

  else if (strcmp(argv[1], "BackwardEuler") == 0) {
    int optn = 0;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &optn) != TCL_OK) {
        opserr << "WARNING integrator BackwardEuler <option> - undefined "
                  "option specified\n";
        return TCL_ERROR;
      }
    }
    theTransientIntegrator = new BackwardEuler(optn);
  }

  else if (strcmp(argv[1], "Newmark") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_Newmark(rt);

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "GimmeMCK") == 0 || strcmp(argv[1], "ZZTop") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GimmeMCK(rt);

  } else if (strcmp(argv[1], "StagedNewmark") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_StagedNewmark(rt);

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);

#ifdef OPS_USE_PFEM
  } else if (strcmp(argv[1], "PFEM") == 0) {
    theTransientIntegrator = new PFEMIntegrator();
#endif

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "NewmarkExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkExplicit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrReduct(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrLimit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "NewmarkHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSFixedNumIter(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHT") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHT_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "GeneralizedAlpha") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GeneralizedAlpha(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "AlphaOS") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "AlphaOS_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized_TP(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "Collocation") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_Collocation(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CollocationHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrReduct(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CollocationHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrLimit(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CollocationHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSFixedNumIter(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "Newmark1") == 0) {
    double gamma;
    double beta;
    double alphaM, betaK, betaKi, betaKc;
    if (argc != 4 && argc != 8) {
      opserr << "WARNING integrator Newmark1 gamma beta <alphaM> "
                "<betaKcurrent> <betaKi> <betaKlastCommitted>\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined gamma\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined beta\n";
      return TCL_ERROR;
    }

    if (argc == 8 || argc == 7) {
      if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - alphaM\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaK\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKi\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKc\n";
        return TCL_ERROR;
      }
    }
    if (argc == 4)
      theTransientIntegrator = new Newmark1(gamma, beta);
    else
      theTransientIntegrator =
          new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc);

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "WilsonTheta") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_WilsonTheta(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "ExplicitDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_ExplicitDifference(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CentralDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifference(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CentralDifferenceAlternative") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceAlternative(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "CentralDifferenceNoDamping") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceNoDamping(rt);

    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "Transient") == 0) {

    theTransientIntegrator = 0;

    // try existing loaded packages
    ExternalClassFunction *integratorCommands =
        theExternalTransientIntegratorCommands;
    bool found = false;
    //    int result = TCL_ERROR;
    while (integratorCommands != NULL && found == false) {

      if (strcmp(argv[2], integratorCommands->funcName) == 0) {

          OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, domain);
        void *theRes = (*(integratorCommands->funcPtr))();
        if (theRes != 0) {
          theTransientIntegrator = (TransientIntegrator *)theRes;
          found = true;
        }
      } else
        integratorCommands = integratorCommands->next;
    }

    //
    // if not there try loading package
    //

    if (found == false) {

      void *libHandle;
      void *(*funcPtr)();
      int integratorNameLength = strlen(argv[2]);
      char *tclFuncName = new char[integratorNameLength + 5];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[2]);

      int res = getLibraryFunction(argv[2], tclFuncName, &libHandle,
                                   (void **)&funcPtr);

      delete[] tclFuncName;

      if (res == 0) {

        char *integratorName = new char[integratorNameLength + 1];
        strcpy(integratorName, argv[2]);
        ExternalClassFunction *theIntegratorCommand = new ExternalClassFunction;
        theIntegratorCommand->funcPtr = funcPtr;
        theIntegratorCommand->funcName = integratorName;
        theIntegratorCommand->next = theExternalTransientIntegratorCommands;
        theExternalTransientIntegratorCommands = theIntegratorCommand;

        OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, domain);

        void *theRes = (*funcPtr)();
        if (theRes != 0) {
          theTransientIntegrator = (TransientIntegrator *)theRes;
        }
      }
    }

    if (theTransientIntegrator == 0) {
      opserr << "Transient Integrator Not Found \n";
      return TCL_ERROR;
    }

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1], "Static") == 0) {

    the_static_integrator = 0;

    // try existing loaded packages
    ExternalClassFunction *integratorCommands =
        theExternalStaticIntegratorCommands;
    bool found = false;

    while (integratorCommands != NULL && found == false) {

      if (strcmp(argv[2], integratorCommands->funcName) == 0) {

          OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, domain);
        void *theRes = (*(integratorCommands->funcPtr))();
        if (theRes != 0) {
          the_static_integrator = (StaticIntegrator *)theRes;
          found = true;
        }
      } else
        integratorCommands = integratorCommands->next;
    }

    //
    // if not there try loading package
    //

    if (found == false) {

      void *libHandle;
      void *(*funcPtr)();
      int integratorNameLength = strlen(argv[2]);
      char *tclFuncName = new char[integratorNameLength + 5];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[2]);

      int res = getLibraryFunction(argv[2], tclFuncName, &libHandle,
                                   (void **)&funcPtr);

      delete[] tclFuncName;

      if (res == 0) {

        char *integratorName = new char[integratorNameLength + 1];
        strcpy(integratorName, argv[2]);
        ExternalClassFunction *theIntegratorCommand = new ExternalClassFunction;
        theIntegratorCommand->funcPtr = funcPtr;
        theIntegratorCommand->funcName = integratorName;
        theIntegratorCommand->next = theExternalStaticIntegratorCommands;
        theExternalStaticIntegratorCommands = theIntegratorCommand;

        OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, domain);

        void *theRes = (*funcPtr)();
        if (theRes != 0) {
          the_static_integrator = (StaticIntegrator *)theRes;
        }
      }
    }

    if (the_static_integrator == 0) {
      opserr << "Static Integrator Not Found \n";
      return TCL_ERROR;
    }

    // if the analysis exists - we want to change the Integrator
    if (the_static_analysis != 0)
      the_static_analysis->setIntegrator(*the_static_integrator);
  }

  else {
    opserr << "WARNING No Integrator type exists \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

