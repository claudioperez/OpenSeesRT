#include <G3Parse.h>
#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>


// AddingSensitivity:BEGIN ////////////////////////////////////
int
sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - sensNodeDisp nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeDisp nodeTag? dof? paramTag?- could not read nodeTag? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeDisp nodeTag? dof? paramTag?- could not read dof? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << "WARNING nodeDisp paramTag? dof? paramTag?- could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  Node *theNode = theDomain->getNode(tag);
  if (theNode == 0) {
    opserr << "sensNodeDisp: node " << tag << " not found" << endln;
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == 0) {
    opserr << "sensNodeDisp: parameter " << paramTag << " not found" << endln;
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getDispSensitivity(dof, gradIndex);

  char buffer[40];
  sprintf(buffer, "%35.20f", value);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - sensNodeVel nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read "
              "dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read "
              "paramTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain->getNode(tag);
  if (theNode == 0) {
    opserr << "sensNodeVel: node " << tag << " not found" << endln;
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == 0) {
    opserr << "sensNodeVel: parameter " << paramTag << " not found" << endln;
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getVelSensitivity(dof, gradIndex);

  char buffer[40];
  sprintf(buffer, "%35.20f", value);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
sensNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - sensNodeAccel nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sensNodeAccel nodeTag? dof? paramTag? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING sensNodeAccel nodeTag? dof? paramTag? - could not read "
              "dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << "WARNING sendNodeAccel nodeTag? dof? paramTag? - could not read "
              "paramTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain->getNode(tag);
  if (theNode == 0) {
    opserr << "sensNodeAccel: node " << tag << " not found" << endln;
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == 0) {
    opserr << "sensNodeAccel: parameter " << paramTag << " not found" << endln;
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getAccSensitivity(dof, gradIndex);

  char buffer[40];
  sprintf(buffer, "%35.20f", value);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
sensNodePressure(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sensNodePressure nodeTag? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sensNodePressure nodeTag? paramTag?- could not read "
              "nodeTag? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << "WARNING sensNodePressure paramTag? paramTag?- could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  double dp = 0.0;
  Pressure_Constraint *thePC = theDomain->getPressure_Constraint(tag);
  if (thePC != 0) {
    // int ptag = thePC->getPressureNode();
    // Node* pNode = theDomain->getNode(ptag);
    Node *pNode = thePC->getPressureNode();
    if (pNode != 0) {

      Parameter *theParam = theDomain->getParameter(paramTag);
      if (theParam == 0) {
        opserr << "sensNodePressure: parameter " << paramTag << " not found"
               << endln;
        return TCL_ERROR;
      }

      int gradIndex = theParam->getGradIndex();
      dp = pNode->getVelSensitivity(1, gradIndex);
    }
  }

  char buffer[40];
  sprintf(buffer, "%35.20f", dp);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
sensSectionForce(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
#ifdef _RELIABILITY
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr
        << "WARNING want - sensSectionForce eleTag? <secNum?> dof? paramTag?\n";
    return TCL_ERROR;
  }

  // opserr << "sensSectionForce: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag, dof, paramTag;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read eleTag? \n";
    return TCL_ERROR;
  }

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 4) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- "
                "could not read secNum? \n";
      return TCL_ERROR;
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &paramTag) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read paramTag? \n";
    return TCL_ERROR;
  }

  ParameterIter &pIter = theDomain->getParameters();
  Parameter *theParam;
  while ((theParam = pIter()) != 0)
    theParam->activate(false);

  theParam = theDomain->getParameter(paramTag);
  int gradIndex = theParam->getGradIndex();
  theParam->activate(true);

  Element *theElement = theDomain->getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sensSectionForce element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "dsdh";
  const char *argvv[3];
  int argcc = 3;
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;
  if (argc < 5) { // For zeroLengthSection
    argcc = 2;
    argvv[1] = c;
  }

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    Tcl_SetResult(interp, "0.0", TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponseSensitivity(gradIndex);
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer, "%12.8g", theVec(dof - 1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  theParam->activate(false);

  delete theResponse;
#endif
  return TCL_OK;
}

///////////////////////Abbas//////////////////////////////

int
sensLambda(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* domain = G3_getDomain(rt);
  StaticAnalysis* the_static_analysis = G3_getStaticAnalysis(rt);
  StaticIntegrator* the_static_integrator = G3_getStaticIntegrator(rt);

  if (argc < 3) {
    opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }

  int pattern, paramTag;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }

  LoadPattern *thePattern = domain->getLoadPattern(pattern);
  if (thePattern == 0) {
    opserr << "ERROR load pattern with tag " << pattern
           << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << "WARNING sensLambda patternTag?  paramTag?- could not read "
              "paramTag? ";
    return TCL_ERROR;
  }
  Parameter *theParam = domain->getParameter(paramTag);
  if (theParam == 0) {
    opserr << "sensLambda: parameter " << paramTag << " not found" << endln;
    return TCL_ERROR;
  }

  IncrementalIntegrator *theIntegrator = 0;

  if (the_static_analysis != 0 && the_static_integrator != 0) {
    theIntegrator = the_static_integrator;
  } else if (theTransientAnalysis != 0 && theTransientIntegrator != 0) {
    theIntegrator = theTransientIntegrator;
  }

  int gradIndex = theParam->getGradIndex();
  double factor = thePattern->getLoadFactorSensitivity(gradIndex);

  char buffer[40];
  sprintf(buffer, "%35.20f", factor);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

////////////////////////////Abbas///////////////////////////////////////



int
setParameter(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theDomain = G3_getDomain(rt);

  int argLoc = 1;
  double newValue = 0.0;
  ID eleIDs(0, 32);
  int numEle = 0;
  int flag = 0;

  if (strstr(argv[argLoc], "-val") != 0) {
    if (Tcl_GetDouble(interp, argv[argLoc + 1], &newValue) != TCL_OK) {
      opserr << "WARNING setParameter: invalid parameter value\n";
      return TCL_ERROR;
    }
  } else {
    opserr << "WARNING setParameter:  -val not found " << endln;
    return TCL_ERROR;
  }

  argLoc += 2;

  if (strstr(argv[argLoc], "-ele") != 0) {

    if ((strcmp(argv[argLoc], "-ele") == 0) ||
        (strcmp(argv[argLoc], "-eles") == 0) ||
        (strcmp(argv[argLoc], "-element") == 0)) {

      //
      // read in a list of ele until end of command or other flag
      //

      argLoc++;
      int eleTag;

      while (argLoc < argc &&
             Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
        eleIDs[numEle] = eleTag;
        numEle++;
        argLoc++;
      }

      if (numEle > 0)
        flag = 1;

    } else if (strcmp(argv[argLoc], "-eleRange") == 0) {

      flag = 2;

      // ensure no segmentation fault if user messes up
      if (argc < argLoc + 3) {
        opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no "
                  "ele tags specified\n";
        return TCL_ERROR;
      }

      //
      // read in start and end tags of two elements & add set [start,end]
      //

      int start, end;
      if (Tcl_GetInt(interp, argv[argLoc + 1], &start) != TCL_OK) {
        opserr
            << "WARNING recorder Element -eleRange start? end? - invalid start "
            << argv[argLoc + 1] << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[argLoc + 2], &end) != TCL_OK) {
        opserr
            << "WARNING recorder Element -eleRange start? end? - invalid end "
            << argv[argLoc + 2] << endln;
        return TCL_ERROR;
      }
      if (start > end) {
        int swap = end;
        end = start;
        start = swap;
      }
      eleIDs[0] = start;
      eleIDs[1] = end;

      argLoc += 3;
    }

    ElementStateParameter theParameter(newValue, &argv[argLoc], argc - argLoc,
                                       flag, &eleIDs);

    theDomain->addParameter(&theParameter);
  }

  return TCL_OK;
}
