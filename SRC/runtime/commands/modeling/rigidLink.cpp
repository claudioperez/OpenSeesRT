/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <assert.h>
#include <tcl.h>

#include <runtime/BasicModelBuilder.h>
#include "RigidRod.h"
#include "RigidBeam.h"
#include "RigidDiaphragm.h"

int
TclCommand_RigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theTclDomain = ((BasicModelBuilder*)clientData)->getDomain();

  if (argc < 3) {
      opserr << "WARNING rigidLink perpDirn? rNode? <cNodes?>\n";
      return TCL_ERROR;
  }

  int rNode, perpDirn;
  if (Tcl_GetInt(interp, argv[1], &perpDirn) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read perpDirn? \n";
      return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read rNode \n";
      return TCL_ERROR;
  }

  // read in the constrained Nodes
  int numConstrainedNodes = argc - 3;
  ID constrainedNodes(numConstrainedNodes);
  for (int i=0; i<numConstrainedNodes; i++) {
      int cNode;
      if (Tcl_GetInt(interp, argv[3+i], &cNode) != TCL_OK) {
          opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read a cNode\n";
          return TCL_ERROR;
      }
      constrainedNodes(i) = cNode;
  }

  RigidDiaphragm theLink(*theTclDomain, rNode, constrainedNodes, perpDirn-1);

  return TCL_OK;
}



int
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theTclDomain = ((BasicModelBuilder*)clientData)->getDomain();

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

