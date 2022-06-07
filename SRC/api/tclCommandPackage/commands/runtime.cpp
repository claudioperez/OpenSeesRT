
#include <g3_api.h>
#include <Domain.h>
#include <tcl.h>

int
TclCommand_setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* domain = G3_getDomain(rt);
  
  domain->setLoadConstant();
  if (argc == 3) {
    if (strcmp(argv[1], "-time") == 0) {
      double newTime;
      if (Tcl_GetDouble(interp, argv[2], &newTime) != TCL_OK) {
        opserr << "WARNING readingvalue - loadConst -time value \n";
        return TCL_ERROR;
      } else {
        domain->setCurrentTime(newTime);
        domain->setCommittedTime(newTime);
      }
    }
  }
  return TCL_OK;
}

int
TclCommand_setCreep(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING illegal command - setCreep value? \n";
    return TCL_ERROR;
  }
  int newFlag;
  if (Tcl_GetInt(interp, argv[1], &newFlag) != TCL_OK) {
    opserr << "WARNING reading creep value - setCreep newFlag? \n";
    return TCL_ERROR;
  } else {
    G3_getDomain(G3_getRuntime(interp))->setCreep(newFlag);
  }
  return TCL_OK;
}

int
TclCommand_setTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* domain = G3_getDomain(rt);
  if (argc < 2) {
    opserr << "WARNING illegal command - time pseudoTime? \n";
    return TCL_ERROR;
  }
  double newTime;
  if (Tcl_GetDouble(interp, argv[1], &newTime) != TCL_OK) {
    opserr << "WARNING reading time value - time pseudoTime? \n";
    return TCL_ERROR;
  } else {
    domain->setCurrentTime(newTime);
    domain->setCommittedTime(newTime);
  }
  return TCL_OK;
}

int
TclCommand_getTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  double time = domain->getCurrentTime();

  // get the display format
  char format[80];
  if (argc == 1) {
    sprintf(format, "%f", time);
  } else if (argc == 2) {
    sprintf(format, argv[1], time);
  }

  // now we copy the value to the tcl string that is returned
  //  sprintf(interp->result,format,time);
  Tcl_SetResult(interp, format, TCL_VOLATILE);
  return TCL_OK;
}

