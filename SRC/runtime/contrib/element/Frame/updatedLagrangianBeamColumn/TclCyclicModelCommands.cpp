//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <CyclicModel.h>
#include <LinearCyclic.h>
#include <BilinearCyclic.h>
#include <QuadraticCyclic.h>
#include <runtime/BasicModelBuilder.h>

int
TclBasicBuilder_addLinearCylic(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << endln;
    return TCL_ERROR;
  }

  CyclicModel *cModel = new LinearCyclic(tag);
  if (!cModel) {
    opserr << "TclBasicBuilder_addLinearCycylic - could not allocate memory\n";
    return TCL_ERROR;
  }
  if (builder->addRegistryObject("CyclicModel", tag, (void*)cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << endln;
    opserr << "\a";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addBilinearCyclic(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;
  double wt;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK) {
    opserr << "WARNING invalid arg[3]" << endln;
    return TCL_ERROR;
  }

  CyclicModel *cModel = new BilinearCyclic(tag, wt);
  if (builder->addRegistryObject("CyclicModel", tag, (void*)cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << endln;
    opserr << "\a";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addQuadraticCyclic(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;
  double wt, qy;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK) {
    opserr << "WARNING invalid arg[3]" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &qy) != TCL_OK) {
    opserr << "WARNING invalid arg[4]" << endln;
    return TCL_ERROR;
  }

  CyclicModel *cModel = new QuadraticCyclic(tag, wt, qy);
  if (builder->addRegistryObject("CyclicModel", tag, (void*)cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << endln;
    opserr << "\a";
    return TCL_ERROR;
  }
  return TCL_OK;
}

/*******************************************************************************************/
int
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "linear") == 0) {
    int result = TclBasicBuilder_addLinearCylic(clientData, interp, argc, argv);
    return result;

  } else if (strcmp(argv[1], "bilinear") == 0) {
    int result = TclBasicBuilder_addBilinearCyclic(clientData, interp, argc, argv);
    return result;
  }

  else if (strcmp(argv[1], "quadratic") == 0) {
    int result = TclBasicBuilder_addQuadraticCyclic(clientData, interp, argc, argv);
    return result;
  }

  else

    return TCL_ERROR;
}
