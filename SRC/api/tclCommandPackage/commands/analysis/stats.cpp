
int totalCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int solveCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int accelCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int numFact(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int numIter(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int systemSize(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);

{
  Tcl_CreateCommand(interp, "totalCPU",   &totalCPU, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "solveCPU",   &solveCPU, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "accelCPU",   &accelCPU, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "numFact",    &numFact, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "numIter",    &numIter, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "systemSize", &systemSize, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
}



int
totalCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getTotalTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
solveCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getSolveTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
accelCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getAccelTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
numFact(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%d", theAlgorithm->getNumFactorizations());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
systemSize(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theSOE == 0) {
    sprintf(buffer, "NO SYSTEM SET");
    return TCL_OK;
  }

  sprintf(buffer, "%d", theSOE->getNumEqn());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
numIter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%d", theAlgorithm->getNumIterations());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

