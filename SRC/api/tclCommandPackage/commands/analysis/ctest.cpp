#include <G3Parse.h>

#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>

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

extern ConvergenceTest *theTest;
extern DirectIntegrationAnalysis *theTransientAnalysis;

ConvergenceTest*
RT_newConvergenceTest(G3_Runtime* rt, int argc, G3_Char** argv);


int
specifyCTest(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << "WARNING need to specify a ConvergenceTest Type type \n";
    return TCL_ERROR;
  }
  G3_Runtime *rt = G3_getRuntime(interp);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);

  ConvergenceTest* theNewTest = RT_newConvergenceTest(rt, argc, argv);
  if (theNewTest != 0) {
    theTest = theNewTest;

    // if the analysis exists - we want to change the Test
    if (the_static_analysis != 0)
      the_static_analysis->setConvergenceTest(*theTest);

    else if (theTransientAnalysis != 0)
      theTransientAnalysis->setConvergenceTest(*theTest);

#ifdef _PARALLEL_PROCESSING
    if (the_static_analysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
        theSub->setAnalysisConvergenceTest(*theTest);
        ;
      }
    }
#endif
  }
}

ConvergenceTest*
RT_newConvergenceTest(G3_Runtime* rt, int argc, G3_Char** argv)
{
  Domain *domain = G3_getDomain(rt);

  // get the tolerence first
  double tol = 0.0;
  double tol2 = 0.0;
  double tolp = 0.0;
  double tolp2 = 0.0;
  double tolrel = 0.0;
  double tolprel = 0.0;
  double maxTol = OPS_MAXTOL;

  int numIter = 0;
  int printIt = 0;
  int normType = 2;
  int maxIncr = -1;

  if ((strcmp(argv[1], "NormDispAndUnbalance") == 0) ||
      (strcmp(argv[1], "NormDispOrUnbalance") == 0)) {
    if (argc == 5) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[3], &tol2) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &numIter) != TCL_OK)
        return nullptr;
    } else if (argc == 6) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[3], &tol2) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[5], &printIt) != TCL_OK)
        return nullptr;
    } else if (argc == 7) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[3], &tol2) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[5], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[6], &normType) != TCL_OK)
        return nullptr;
    } else if (argc == 8) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[3], &tol2) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[5], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[6], &normType) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[7], &maxIncr) != TCL_OK)
        return nullptr;
    }

#ifdef OPS_USE_PFEM
  } else if (strcmp(argv[1], "PFEM") == 0) {
    if (argc > 8) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[3], &tolp) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[4], &tol2) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[5], &tolp2) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[6], &tolrel) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[7], &tolprel) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[8], &numIter) != TCL_OK)
        return nullptr;
    }
    if (argc > 9) {
      if (G3Parse_getInt(rt, argv[9], &maxIncr) != TCL_OK)
        return nullptr;
    }
    if (argc > 10) {
      if (G3Parse_getInt(rt, argv[10], &printIt) != TCL_OK)
        return nullptr;
    }
    if (argc > 11) {
      if (G3Parse_getInt(rt, argv[11], &normType) != TCL_OK)
        return nullptr;
    }
#endif
  } else if (strcmp(argv[1], "FixedNumIter") == 0) {

    if (argc == 3) {
      if (G3Parse_getInt(rt, argv[2], &numIter) != TCL_OK)
        return nullptr;
    } else if (argc == 4) {
      if (G3Parse_getInt(rt, argv[2], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &printIt) != TCL_OK)
        return nullptr;
    } else if (argc == 5) {
      if (G3Parse_getInt(rt, argv[2], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &normType) != TCL_OK)
        return nullptr;
    } else if (argc == 6) {
      if (G3Parse_getInt(rt, argv[2], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &normType) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[5], &maxTol) != TCL_OK)
        return nullptr;
    }

  } else {
    if (argc == 4) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
    } else if (argc == 5) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &printIt) != TCL_OK)
        return nullptr;
    } else if (argc == 6) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[5], &normType) != TCL_OK)
        return nullptr;
    } else if (argc == 7) {
      if (G3Parse_getDouble(rt, argv[2], &tol) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[4], &printIt) != TCL_OK)
        return nullptr;
      if (G3Parse_getInt(rt, argv[5], &normType) != TCL_OK)
        return nullptr;
      if (G3Parse_getDouble(rt, argv[6], &maxTol) != TCL_OK)
        return nullptr;
    }
  }

  ConvergenceTest *theNewTest = 0;

  if (numIter == 0) {
    opserr << "ERROR: no numIter specified in test command\n";
    return nullptr;
  }

  if (strcmp(argv[1], "FixedNumIter") == 0)
    theNewTest = new CTestFixedNumIter(numIter, printIt, normType);
  else {
    if (tol == 0.0) {
      opserr << "ERROR: no tolerance specified in test command\n";
      return nullptr;
    }
    if (strcmp(argv[1], "NormUnbalance") == 0)
      theNewTest = new CTestNormUnbalance(tol, numIter, printIt, normType,
                                          maxIncr, maxTol);
    else if (strcmp(argv[1], "NormDispIncr") == 0)
      theNewTest =
          new CTestNormDispIncr(tol, numIter, printIt, normType, maxTol);
    else if (strcmp(argv[1], "NormDispAndUnbalance") == 0)
      theNewTest = new NormDispAndUnbalance(tol, tol2, numIter, printIt,
                                            normType, maxIncr);
    else if (strcmp(argv[1], "NormDispOrUnbalance") == 0)
      theNewTest = new NormDispOrUnbalance(tol, tol2, numIter, printIt,
                                           normType, maxIncr);
    else if (strcmp(argv[1], "EnergyIncr") == 0)
      theNewTest = new CTestEnergyIncr(tol, numIter, printIt, normType, maxTol);
    else if (strcmp(argv[1], "RelativeNormUnbalance") == 0)
      theNewTest =
          new CTestRelativeNormUnbalance(tol, numIter, printIt, normType);
    else if (strcmp(argv[1], "RelativeNormDispIncr") == 0)
      theNewTest =
          new CTestRelativeNormDispIncr(tol, numIter, printIt, normType);
    else if (strcmp(argv[1], "RelativeEnergyIncr") == 0)
      theNewTest = new CTestRelativeEnergyIncr(tol, numIter, printIt, normType);
    else if (strcmp(argv[1], "RelativeTotalNormDispIncr") == 0)
      theNewTest =
          new CTestRelativeTotalNormDispIncr(tol, numIter, printIt, normType);
#ifdef OPS_USE_PFEM
    else if (strcmp(argv[1], "PFEM") == 0)
      theNewTest = new CTestPFEM(tol, tolp, tol2, tolp2, tolrel, tolprel,
                                 numIter, maxIncr, printIt, normType);
#endif
    else {
      opserr << "WARNING No ConvergenceTest type (NormUnbalance, NormDispIncr, "
                "EnergyIncr, \n";
      opserr << "RelativeNormUnbalance, RelativeNormDispIncr, "
                "RelativeEnergyIncr, \n";
      opserr << "RelativeTotalNormDispIncr, FixedNumIter)\n";
      return nullptr;
    }
  }

  return theNewTest;
}

int
getCTestNorms(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  if (theTest != 0) {
    const Vector &data = theTest->getNorms();

    char buffer[40];
    int size = data.Size();
    for (int i = 0; i < size; i++) {
      sprintf(buffer, "%35.20e", data(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }

    return TCL_OK;
  }

  opserr << "ERROR testNorms - no convergence test!\n";
  return TCL_ERROR;
}

int
getCTestIter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTest != 0) {
    int res = theTest->getNumTests();

    char buffer[10];
    sprintf(buffer, "%d", res);
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
  }

  opserr << "ERROR testIter - no convergence test!\n";
  return TCL_ERROR;
}

