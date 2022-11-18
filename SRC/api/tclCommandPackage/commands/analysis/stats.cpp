
// int totalCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int solveCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int accelCPU(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
int numFact(ClientData  clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
// int numIter(ClientData  clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);
// int systemSize(ClientData clientData, Tcl_Interp *interp, int argc,TCL_Char **argv);

{
//Tcl_CreateCommand(interp, "totalCPU",   &totalCPU,   (ClientData)NULL, nullptr);
//Tcl_CreateCommand(interp, "solveCPU",   &solveCPU,   (ClientData)NULL, nullptr);
  Tcl_CreateCommand(interp, "accelCPU",   &accelCPU,   (ClientData)NULL, nullptr);
  Tcl_CreateCommand(interp, "numFact",    &numFact,    (ClientData)NULL, nullptr);
//Tcl_CreateCommand(interp, "numIter",    &numIter,    (ClientData)NULL, nullptr);
//Tcl_CreateCommand(interp, "systemSize", &systemSize, (ClientData)NULL, nullptr);
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


