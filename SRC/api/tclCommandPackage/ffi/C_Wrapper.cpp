#include <g3_api.h>
#include <tcl.h>
#include <Domain.h>
#include <UniaxialMaterial.h>
#include <WrapperUniaxialMaterial.h>

// extern TCL_Char **currentArgv;
// extern int currentArg;
// extern int maxArg;

UniaxialMaterial *
Tcl_addWrapperUniaxialMaterial(matObj *theMat, ClientData clientData,
                               Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // cmp //  theInterp = interp;

  // currentArgv = argv;
  // currentArg = 2;
  // maxArg = argc;

  G3_Runtime* rt = G3_getRuntime(interp);
  Domain* theDomain = G3_getDomain(rt);

  // get the current load factor
  static modelState theModelState;
  if (theDomain != 0) {
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;
    theModelState.time = time;
    theModelState.dt = dt;
  }

  // invoke the mat function with isw = 0
  int isw = ISW_INIT;
  int result = 0;
  theMat->matFunctPtr(theMat, &theModelState, 0, 0, 0, &isw, &result);
  int matType = theMat->matType; // GR added to support material

  if (result != 0 || matType != OPS_UNIAXIAL_MATERIAL_TYPE) {
    opserr << "Tcl_addWrapperUniaxialMaterial - failed in element function "
           << result << endln;
    return NULL;
  }

  WrapperUniaxialMaterial *theMaterial =
      new WrapperUniaxialMaterial(argv[1], theMat);

  return theMaterial;
}

/*
extern "C" int
OPS_AllocateMaterial(matObject *theMat)
{

  if (theMat->nParam > 0)
    theMat->theParam = new double[theMat->nParam];

  int nState = theMat->nState;

  if (nState > 0) {
    theMat->cState = new double[nState];
    theMat->tState = new double[nState];
    for (int i = 0; i < nState; i++) {
      theMat->cState[i] = 0;
      theMat->tState[i] = 0;
    }
  } else {
    theMat->cState = 0;
    theMat->tState = 0;
  }

  return 0;
}

extern "C" int
OPS_AllocateElement(eleObject *theEle, int *matTags, int *matType)
{
  if (theEle->nNode > 0)
    theEle->node = new int[theEle->nNode];

  if (theEle->nParam > 0)
    theEle->param = new double[theEle->nParam];

  if (theEle->nState > 0) {
    theEle->cState = new double[theEle->nState];
    theEle->tState = new double[theEle->nState];
  }

  int numMat = theEle->nMat;
  if (numMat > 0)
    theEle->mats = new matObject *[numMat];

  for (int i = 0; i < numMat; i++) {

    matObject *theMat = OPS_GetMaterial(&(matTags[i]), matType);
    //    matObject *theMat = OPS_GetMaterial(&(matTags[i]));

    theEle->mats[i] = theMat;
  }

  return 0;
}


extern "C" int
OPS_GetNodeCrd(int *nodeTag, int *sizeCrd, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);
  if (theNode == 0) {
    opserr << "OPS_GetNodeCrd - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeCrd;
  const Vector &crd = theNode->getCrds();
  if (crd.Size() != size) {
    opserr << "OPS_GetNodeCrd - crd size mismatch\n";
    opserr << "Actual crd size is: " << crd.Size()
           << endln; // MRL Add Error Detection
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = crd(i);

  return 0;
}

extern "C" int
OPS_GetNodeDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeDisp - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getTrialDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeDisp - crd size mismatch\n";
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = disp(i);

  return 0;
}

extern "C" int
OPS_GetNodeVel(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeVel - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &vel = theNode->getTrialVel();

  if (vel.Size() != size) {
    opserr << "OPS_GetNodeVel - crd size mismatch\n";
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = vel(i);

  return 0;
}

extern "C" int
OPS_GetNodeAccel(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeAccel - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &accel = theNode->getTrialAccel();

  if (accel.Size() != size) {
    opserr << "OPS_GetNodeAccel - accel size mismatch\n";
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = accel(i);

  return 0;
}

extern "C" int
OPS_GetNodeIncrDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeIncrDisp - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getIncrDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeIncrDis - crd size mismatch\n";
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = disp(i);

  return 0;
}

extern "C" int
OPS_GetNodeIncrDeltaDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeIncrDeltaDisp - no node with tag " << *nodeTag
           << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getIncrDeltaDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeIncrDis - crd size mismatch\n";
    return -1;
  }
  for (int i = 0; i < size; i++)
    data[i] = disp(i);

  return 0;
}
*/

int
Tcl_addWrapperElement(eleObj *theEle, ClientData clientData, Tcl_Interp *interp,
                      int argc, TCL_Char **argv, Domain *theDomain,
                      TclBuilder *theModelBuilder)
{
  theInterp = interp;
  // theDomain = domain;
  // theModelBuilder = builder;
  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

  // get the current load factor
  double time = theDomain->getCurrentTime();
  double dt = theDomain->getCurrentTime() - time;

  static modelState theModelState;
  theModelState.time = time;
  theModelState.dt = dt;

  // invoke the ele function with isw = 0
  int isw = ISW_INIT;
  int result = 0;
  theEle->eleFunctPtr(theEle, &theModelState, 0, 0, &isw, &result);

  if (result != 0) {
    opserr << "Tcl_addWrapperElement - failed in element function " << result
           << endln;
    return TCL_ERROR;
  }

  WrapperElement *theElement = new WrapperElement(argv[1], theEle);

  if (theDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element of type: " << argv[1]
           << " to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  return 0;
}




extern "C" int
OPS_InvokeMaterial(eleObject *theEle, int *mat, modelState *model,
                   double *strain, double *stress, double *tang, int *isw)
{
  int error = 0;

  matObject *theMat = theEle->mats[*mat];
  /* fprintf(stderr,"invokeMaterial Address %d %d %d\n",*mat, theMat,
   * sizeof(int)); */

  if (theMat != 0)
    theMat->matFunctPtr(theMat, model, strain, tang, stress, isw, &error);
  else
    error = -1;

  return error;
}

extern "C" int
OPS_InvokeMaterialDirectly(matObject **theMat, modelState *model,
                           double *strain, double *stress, double *tang,
                           int *isw)
{
  int error = 0;
  //  fprintf(stderr,"invokeMaterialDirectly Address %d %d %d\n",theMat,
  //  sizeof(int), *theMat);
  if (*theMat != 0)
    (*theMat)->matFunctPtr(*theMat, model, strain, tang, stress, isw, &error);
  else
    error = -1;

  return error;
}

