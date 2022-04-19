
#include <tcl.h>
#include <RigidRod.h>
#include <RigidBeam.h>
#include <elementAPI.h>

int
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theTclDomain = G3_getDomain(rt);

  if (argc < 4) {
      opserr << "WARNING rigidLink linkType? rNode? cNode?\n";
      return TCL_ERROR;
  }

  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read rNode \n"; 
      return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read CNode \n"; 
      return TCL_ERROR;
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1],"-bar") == 0) || (strcmp(argv[1],"bar") == 0)) {
    RigidRod theLink(*theTclDomain, rNode, cNode);
  } else if ((strcmp(argv[1],"-beam") == 0) || (strcmp(argv[1],"beam") == 0)) {
    RigidBeam theLink(*theTclDomain, rNode, cNode);
  } else {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) \n"; 
      return TCL_ERROR;
  }
  return TCL_OK;
}

