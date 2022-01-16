#include <math.h>

#include <g3_api.h>
#include <tcl.h>
#include <OPS_Globals.h>

#include <Vector.h>
#include <Domain.h>
#include <Node.h>

# include <StandardStream.h>
# include <FileStream.h>
# include <DummyStream.h>
#include <iostream>

int
setNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - setNodeAccel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING setNodeAccel -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read "
              "value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getAccel();
    vel(dof) = value;
    theNode->setTrialAccel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

int
nodeAccel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeAccel nodeTag? dof?\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeAccel nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeAccel nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = the_domain->getNodeResponse(tag, Accel);
  if (nodalResponse == 0)
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {
    if (size < dof)
      return TCL_ERROR;

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned
    // sprintf(interp->result,"%35.20f",value);
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  } else {
    char buffer[40];
    for (int i = 0; i < size; i++) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
nodeResponse(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
    return TCL_ERROR;
  }

  int tag, dof, responseID;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? - could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &responseID) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? responseID? - could not read "
              "responseID? \n";
    return TCL_ERROR;
  }

  dof--;

  const Vector *nodalResponse =
      the_domain->getNodeResponse(tag, (NodeResponseType)responseID);
  if (nodalResponse == 0 || nodalResponse->Size() < dof || dof < 0)
    return TCL_ERROR;

  double value = (*nodalResponse)(dof);

  // now we copy the value to the tcl string that is returned
  //    sprintf(interp->result,"%35.20f",value);
  char buffer[40];
  sprintf(buffer, "%35.20f", value);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
calculateNodalReactions(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain type of system
  int incInertia = 0;

  if (argc == 2) {
    if ((strcmp(argv[1], "-incInertia") == 0) ||
        (strcmp(argv[1], "-dynamical") == 0) ||
        (strcmp(argv[1], "-Dynamic") == 0) ||
        (strcmp(argv[1], "-dynamic") == 0))

      incInertia = 1;

    else if ((strcmp(argv[1], "-rayleigh") == 0))

      incInertia = 2;
  }

  the_domain->calculateNodalReactions(incInertia);

  return TCL_OK;
}

