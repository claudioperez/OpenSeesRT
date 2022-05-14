#include <g3_api.h>
#include <G3Parse.h>
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
#include <ArcLength1.h>
#include <DisplacementControl.h>

#include <Newmark.h>
#include <StagedNewmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
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

StaticIntegrator*
G3Parse_newHSIntegrator(G3_Builder *, int, const char **);
StaticIntegrator*
G3Parse_newLoadControl(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newEQPathIntegrator(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newArcLengthIntegrator(G3_Builder *, int argc, const char *argv[]);
StaticIntegrator*
G3Parse_newStagedLoadControlIntegrator(G3_Builder*, int, TCL_Char **);
StaticIntegrator*
G3Parse_newMinUnbalDispNormIntegrator(G3_Runtime*, int, G3_Char **);
StaticIntegrator*
G3Parse_newDisplacementControlIntegrator(G3_Builder *, int, G3_Char**);
StaticIntegrator*
G3Parse_newStaticIntegrator(G3_Builder *, int, TCL_Char **);

TransientIntegrator*
G3Parse_newNewmark1Integrator(G3_Builder *, int, TCL_Char **);
TransientIntegrator*
G3Parse_newNewmarkIntegrator(G3_Runtime*, int, G3_Char**);
TransientIntegrator*
G3Parse_newTransientIntegrator(G3_Builder *, int, TCL_Char **);
//
// command invoked to allow the Integrator object to be built
//

int
specifyIntegrator(ClientData cd, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  bool assign_to_static_analysis = false;
  G3_Builder *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);
  StaticIntegrator* the_static_integrator = G3_getStaticIntegrator(rt);
  OPS_ResetInputNoBuilder(cd, interp, 2, argc, argv, domain);

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify an Integrator type \n";
    return TCL_ERROR;
  }

  the_static_integrator = G3Parse_newStaticIntegrator(rt, argc, argv);
  if (the_static_analysis != 0) {
    G3_setStaticIntegrator(rt, the_static_integrator);
    the_static_analysis->setIntegrator(*the_static_integrator);
  }

  theTransientIntegrator = G3Parse_newTransientIntegrator(rt, argc, argv);
  if (theTransientAnalysis != 0)
    theTransientAnalysis->setIntegrator(*theTransientIntegrator);

  return TCL_OK;
}

StaticIntegrator*
G3Parse_newStaticIntegrator(G3_Builder* rt, int argc, TCL_Char **argv)
{

  StaticIntegrator* the_static_integrator = nullptr;

  // check argv[1] for type of Numberer and create the object

  if (strcmp(argv[1], "LoadControl") == 0) {
    the_static_integrator = G3Parse_newLoadControl(rt, argc, argv);

  } else if (strcmp(argv[1], "StagedLoadControl") == 0) {
    the_static_integrator = G3Parse_newStagedLoadControlIntegrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "EQPath") == 0) {
    the_static_integrator = G3Parse_newEQPathIntegrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "ArcLength") == 0) {
    the_static_integrator = G3Parse_newArcLengthIntegrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "MinUnbalDispNorm") == 0) {
    the_static_integrator = G3Parse_newMinUnbalDispNormIntegrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "DisplacementControl") == 0) {
    the_static_integrator = G3Parse_newDisplacementControlIntegrator(rt, argc, argv);
  }


  else if (strcmp(argv[1], "ArcLength1") == 0) {
    double arcLength;
    double alpha;
    if (argc != 4) {
      opserr << "WARNING integrator ArcLength1 arcLength alpha \n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[2], &arcLength) != TCL_OK)
      return nullptr;
    if (G3Parse_getDouble(rt, argv[3], &alpha) != TCL_OK)
      return nullptr;
    the_static_integrator = new ArcLength1(arcLength, alpha);
  }

/*
  else if (strcmp(argv[1], "HSConstraint") == 0) {
    the_static_integrator = G3Parse_newHSIntegrator(rt, argc, argv);
  }
*/


  return the_static_integrator;
}

TransientIntegrator*
G3Parse_newTransientIntegrator(G3_Builder *rt, int argc, TCL_Char **argv)
{
  printCommand(argc, argv);
  // G3_Builder *rt = G3_getRuntime(interp);
  if ((strcmp(argv[1], "TRBDF2") == 0) ||
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
      if (G3Parse_getInt(rt, argv[2], &optn) != TCL_OK) {
        opserr << "WARNING integrator BackwardEuler <option> - undefined "
                  "option specified\n";
        return nullptr;
      }
    }
    theTransientIntegrator = new BackwardEuler(optn);
  }

  else if (strcmp(argv[1], "Newmark") == 0) {
    // theTransientIntegrator = (TransientIntegrator *)OPS_Newmark(rt);
    theTransientIntegrator = (TransientIntegrator *)G3Parse_newNewmarkIntegrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "GimmeMCK") == 0 || strcmp(argv[1], "ZZTop") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GimmeMCK(rt);
  } 

  else if (strcmp(argv[1], "StagedNewmark") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_StagedNewmark(rt);
  }

  else if (strcmp(argv[1], "NewmarkExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkExplicit(rt);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrReduct(rt);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrLimit(rt);
  }

  else if (strcmp(argv[1], "NewmarkHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSFixedNumIter(rt);
  }

  else if (strcmp(argv[1], "HHT") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT(rt);
  }

  else if (strcmp(argv[1], "HHT_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT_TP(rt);
  }

  else if (strcmp(argv[1], "HHTGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized(rt);
  }

  else if (strcmp(argv[1], "HHTGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized_TP(rt);
  }

  else if (strcmp(argv[1], "HHTExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit(rt);
  }

  else if (strcmp(argv[1], "HHTExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit_TP(rt);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit(rt);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit_TP(rt);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit(rt);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit_TP(rt);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct(rt);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct_TP(rt);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter(rt);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter_TP(rt);
  }

  else if (strcmp(argv[1], "GeneralizedAlpha") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GeneralizedAlpha(rt);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit(rt);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit_TP(rt);
  }

  else if (strcmp(argv[1], "AlphaOS") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS(rt);
  }

  else if (strcmp(argv[1], "AlphaOS_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS_TP(rt);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized(rt);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized_TP(rt);
  }

  else if (strcmp(argv[1], "Collocation") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_Collocation(rt);
  }

  else if (strcmp(argv[1], "CollocationHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrReduct(rt);
  }

  else if (strcmp(argv[1], "CollocationHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrLimit(rt);
  }

  else if (strcmp(argv[1], "CollocationHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSFixedNumIter(rt);
  }

  else if (strcmp(argv[1], "Newmark1") == 0) {
    theTransientIntegrator = G3Parse_newNewmark1Integrator(rt, argc, argv);
  }

  else if (strcmp(argv[1], "WilsonTheta") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_WilsonTheta(rt);
  }

  else if (strcmp(argv[1], "ExplicitDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_ExplicitDifference(rt);
  }

  else if (strcmp(argv[1], "CentralDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifference(rt);
  }

  else if (strcmp(argv[1], "CentralDifferenceAlternative") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceAlternative(rt);
  }

  else if (strcmp(argv[1], "CentralDifferenceNoDamping") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceNoDamping(rt);

  }

  return theTransientIntegrator;
}

#include <HSConstraint.h>
StaticIntegrator*
G3Parse_newHSIntegrator(G3_Builder* rt, int argc, const char *argv[])
{
    // Tcl_Interp  *interp = G3_getInterpreter(rt);
    double arcLength, psi_u, psi_f, u_ref;

    if (argc < 3) {
      opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> "
                "<u_ref> \n";
      return nullptr;
    }
    if (argc >= 3 && G3Parse_getDouble(rt, argv[2], &arcLength) != TCL_OK)
      return nullptr;
    if (argc >= 4 && G3Parse_getDouble(rt, argv[3], &psi_u) != TCL_OK)
      return nullptr;
    if (argc >= 5 && G3Parse_getDouble(rt, argv[4], &psi_f) != TCL_OK)
      return nullptr;
    if (argc == 6 && G3Parse_getDouble(rt, argv[5], &u_ref) != TCL_OK)
      return nullptr;

    switch (argc) {
    case 3:
      return new HSConstraint(arcLength);
    case 4:
      return new HSConstraint(arcLength, psi_u);
    case 5:
      return new HSConstraint(arcLength, psi_u, psi_f);
    case 6:
      return new HSConstraint(arcLength, psi_u, psi_f, u_ref);
    }
}

StaticIntegrator*
G3Parse_newLoadControl(G3_Builder* rt, int argc, const char *argv[])
{
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator LoadControl dlam <Jd "
                "dlamMin dlamMax>\n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[2], &dLambda) != TCL_OK)
      return nullptr;
    if (argc > 5) {
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[4], &minIncr) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[5], &maxIncr) != TCL_OK)
        return nullptr;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    return new LoadControl(dLambda, numIter, minIncr, maxIncr);
}

#include <EQPath.h>
StaticIntegrator*
G3Parse_newEQPathIntegrator(G3_Builder* rt, int argc, const char *argv[])
{
  // else if (strcmp(argv[1], "EQPath") == 0) {
    // Tcl_Interp  *interp = G3_getInterpreter(rt);
    double arcLength;
    int type;
    int numIter;
    if (argc != 4) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << "REFS : \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
      return nullptr;
    }

    if (G3Parse_getDouble(rt, argv[2], &arcLength) != TCL_OK) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
      return nullptr;
    }

    if (G3Parse_getInt(rt, argv[3], &type) != TCL_OK) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << "$type = 1 Minimum Residual Displacement \n";
      opserr << "$type = 2 Normal Plain \n";
      opserr << "$type = 3 Update Normal Plain \n";
      opserr << "$type = 4 Cylindrical Arc-Length \n";
      return nullptr;
    }

    return new EQPath(arcLength, type);
}

#include <ArcLength.h>
StaticIntegrator *
G3Parse_newArcLengthIntegrator(G3_Builder* rt, int argc, TCL_Char **argv)
{
  // Tcl_Interp *interp = G3_getInterpreter(rt);
  double arcLength;
  double alpha;
  if (argc != 4) {
    opserr << "WARNING integrator ArcLength arcLength alpha \n";
    return nullptr;
  }
  if (G3Parse_getDouble(rt, argv[2], &arcLength) != TCL_OK)
    return nullptr;

  if (G3Parse_getDouble(rt, argv[3], &alpha) != TCL_OK)
    return nullptr;

  return new ArcLength(arcLength, alpha);
}

#include <MinUnbalDispNorm.h>
StaticIntegrator*
G3Parse_newMinUnbalDispNormIntegrator(G3_Runtime* rt, int argc, G3_Char **argv)
{
    double lambda11, minlambda, maxlambda;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j "
                "maxLambda1j>\n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[2], &lambda11) != TCL_OK)
      return nullptr;
    if (argc > 5) {
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[4], &minlambda) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[5], &maxlambda) != TCL_OK)
        return nullptr;
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

    return new MinUnbalDispNorm(lambda11, numIter, minlambda,
                                               maxlambda, signFirstStepMethod);
}

StaticIntegrator*
G3Parse_newDisplacementControlIntegrator(G3_Builder *rt, int argc, G3_Char** argv)
{
    // G3_Builder *rt = G3_getRuntime(interp);
    Domain *domain = G3_getDomain(rt);
    int node, dof, numIter;
    double increment, minIncr, maxIncr;

    if (argc < 5) {
      opserr << "WARNING integrator DisplacementControl node dof dU \n";
      opserr << "<Jd minIncrement maxIncrement>\n";
      return nullptr;
    }
    int tangFlag = 0;

    if (G3Parse_getInt(rt, argv[2], &node) != TCL_OK)
      return nullptr;
    if (G3Parse_getInt(rt, argv[3], &dof) != TCL_OK)
      return nullptr;
    if (G3Parse_getDouble(rt, argv[4], &increment) != TCL_OK)
      return nullptr;

    if (argc == 6 || argc == 9)
      if (argc == 6) {
        if (strcmp(argv[5], "-initial") == 0)
          tangFlag = 1;
      } else if (strcmp(argv[8], "-initial") == 0)
        tangFlag = 1;

    if (argc > 6) {
      if (G3Parse_getInt(rt, argv[5], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[6], &minIncr) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[7], &maxIncr) != TCL_OK)
        return nullptr;
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
      return nullptr;
    }

    int numDOF = theNode->getNumberDOF();
    if (dof <= 0 || dof > numDOF) {
      opserr << "WARNING integrator DisplacementControl node dof dU : invalid "
                "dof given\n";
      return nullptr;
    }

    return
        new DisplacementControl(node, dof-1, increment, domain, numIter,
                                minIncr, maxIncr, tangFlag);

#endif
}

StaticIntegrator*
G3Parse_newStagedLoadControlIntegrator(G3_Builder* rt, int argc, TCL_Char **argv)
{
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator StagedLoadControl dlam "
                "<Jd dlamMin dlamMax>\n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[2], &dLambda) != TCL_OK)
      return nullptr;
    if (argc > 5) {
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[4], &minIncr) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[5], &maxIncr) != TCL_OK)
        return nullptr;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    return  new StagedLoadControl(dLambda, numIter, minIncr, maxIncr);
}

#include <Newmark1.h>
TransientIntegrator*
G3Parse_newNewmark1Integrator(G3_Builder* rt, int argc, TCL_Char **argv)
{
    double gamma;
    double beta;
    double alphaM, betaK, betaKi, betaKc;
    if (argc != 4 && argc != 8) {
      opserr << "WARNING integrator Newmark1 gamma beta <alphaM> "
                "<betaKcurrent> <betaKi> <betaKlastCommitted>\n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[2], &gamma) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined gamma\n";
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[3], &beta) != TCL_OK) {
      opserr << "WARNING integrator Newmark1 gamma beta - undefined beta\n";
      return nullptr;
    }

    if (argc == 8 || argc == 7) {
      if (G3Parse_getDouble(rt, argv[4], &alphaM) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - alphaM\n";
        return nullptr;
      }
      if (G3Parse_getDouble(rt, argv[5], &betaK) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaK\n";
        return nullptr;
      }
      if (G3Parse_getDouble(rt, argv[6], &betaKi) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKi\n";
        return nullptr;
      }
      if (G3Parse_getDouble(rt, argv[7], &betaKc) != TCL_OK) {
        opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi "
                  "betaKc - betaKc\n";
        return nullptr;
      }
    }
    if (argc == 4)
      return new Newmark1(gamma, beta);
    else
      return new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc);
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
      return nullptr;
    }
    if (G3Parse_getInt(rt, argv[2], &node) != TCL_OK)
      return nullptr;
    if (G3Parse_getInt(rt, argv[3], &dof) != TCL_OK)
      return nullptr;
    if (G3Parse_getDouble(rt, argv[4], &increment) != TCL_OK)
      return nullptr;
    if (argc > 7) {
      if (G3Parse_getInt(rt, argv[5], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[6], &minIncr) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[7], &maxIncr) != TCL_OK)
        return nullptr;
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
  }
#endif
