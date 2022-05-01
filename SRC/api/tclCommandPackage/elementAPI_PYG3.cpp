/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
** ****************************************************************** */

/*
** Written: cmp
**/

#include <elementAPI.h>
#include <stdlib.h>
#include <packages.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <g3_api.h>
#include <G3_Runtime.h>
#include <G3_Logging.h>
extern OPS_Stream* opswrnPtr;

#include <TclBasicBuilder.h>
#include <TclSafeBuilder.h>
#include <WrapperElement.h>

#include <map>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <FrictionModel.h>
#include <WrapperUniaxialMaterial.h>
#include <WrapperNDMaterial.h>
#include <LimitCurve.h>
#include <WrapperLimitCurve.h>

#include <DirectIntegrationAnalysis.h>
#include <StaticAnalysis.h>

// use ProfileSPD for default SOE constructor
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>

#include <OPS_Globals.h>

typedef struct elementFunction {
  char *funcName;
  eleFunct theFunct;
  struct elementFunction *next;
} ElementFunction;

typedef struct materialFunction {
  char *funcName;
  matFunct theFunct;
  struct materialFunction *next;
} MaterialFunction;

typedef struct limitCurveFunction {
  char *funcName;
  limCrvFunct theFunct;
  struct limitCurveFunction *next;
} LimitCurveFunction;

extern AnalysisModel *theAnalysisModel;
extern EquiSolnAlgo *theAlgorithm;
extern ConstraintHandler *theHandler;
extern DOF_Numberer *theNumberer;
extern LinearSOE *theSOE;
extern EigenSOE *theEigenSOE;
extern StaticAnalysis *theStaticAnalysis;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern VariableTimeStepDirectIntegrationAnalysis
    *theVariableTimeStepTransientAnalysis;
extern int numEigen;
extern StaticIntegrator *theStaticIntegrator;
extern TransientIntegrator *theTransientIntegrator;
extern ConvergenceTest *theTest;
extern bool builtModel;

static ElementFunction *theElementFunctions = NULL;
static MaterialFunction *theMaterialFunctions = NULL;
static LimitCurveFunction *theLimitCurveFunctions = NULL;

static Tcl_Interp *theInterp = 0;
static Domain *theDomain = 0;

static TclBuilder *theModelBuilder = 0;

static TCL_Char **currentArgv = 0;
static int currentArg = 0;
static int maxArg = 0;

extern const char *getInterpPWD(Tcl_Interp *interp);
extern FE_Datastore *theDatabase;

// static int uniaxialMaterialObjectCount = 0;

modelState theModelState;

struct cmp_str {
  bool
  operator()(const char *a, const char *b)
  {
    return strcmp(a, b) < 0;
  }
};

// std::map<char *, eleFunct, cmp_str>
//     theEleFunctions; // map of user added ele functions
// std::map<char *, eleFunct, cmp_str>
//     theUniaxialMaterialFunctions; // map of user added material functions

// std::map<int, UniaxialMaterial *>theUniaxialMaterials;           // map for
// UniaxialMaterial objects needed by user added ele functions'

static void
OPS_InvokeMaterialObject(struct matObject *theMat, modelState *theModel,
                         double *strain, double *tang, double *stress, int *isw,
                         int *result)
{
  int matType = (int)theMat->theParam[0];

  if (matType == 1) {
    //  UniaxialMaterial *theMaterial = theUniaxialMaterials[matCount];
    UniaxialMaterial *theMaterial = (UniaxialMaterial *)theMat->matObjectPtr;
    if (theMaterial == 0) {
      *result = -1;
      return;
    }

    if (*isw == ISW_COMMIT) {
      *result = theMaterial->commitState();
      return;
    } else if (*isw == ISW_REVERT) {
      *result = theMaterial->revertToLastCommit();
      return;
    } else if (*isw == ISW_REVERT_TO_START) {
      *result = theMaterial->revertToStart();
      return;
    } else if (*isw == ISW_FORM_TANG_AND_RESID) {
      double matStress = 0.0;
      double matTangent = 0.0;
      int res = theMaterial->setTrial(strain[0], matStress, matTangent);
      stress[0] = matStress;
      tang[0] = matTangent;
      *result = res;
      return;
    }
  }

  return;
}

extern "C" int
OPS_Error(char *errorMessage, int length)
{
  opserr << errorMessage;
  opserr << endln;

  return 0;
}

extern "C" int
OPS_GetNumRemainingInputArgs()
{
  return maxArg - currentArg;
}

extern "C" int
OPS_ResetCurrentInputArg(int cArg)
{
  if (cArg < 0)
    currentArg += cArg;
  else
    currentArg = cArg;

  return 0;
}


// extern "C"
int
OPS_ResetInput(ClientData clientData, Tcl_Interp *interp, int cArg, int mArg,
               TCL_Char **argv, Domain *domain, TclBuilder *builder)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  G3_setDomain(rt, domain);
  G3_setModelBuilder(rt, builder);
  theInterp = interp;
  // theDomain = domain;
  theModelBuilder = builder;
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;

  return 0;
}


extern "C" int
OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                        int mArg, TCL_Char **argv, Domain *domain)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  G3_setDomain(rt, domain);
  theInterp = interp;
  // theDomain = domain;
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;
  return 0;
}

extern "C" int
OPS_GetIntInput(int *numData, int *data)
{
  int size = *numData;

  for (int i = 0; i < size; i++) {
    if ((currentArg >= maxArg) ||
        (Tcl_GetInt(theInterp, currentArgv[currentArg], &data[i]) != TCL_OK)) {
      // opserr << "OPS_GetIntInput -- error reading " << currentArg << endln;
      return -1;
    } else
      currentArg++;
  }

  return 0;
}

extern "C" int
OPS_SetIntOutput(int *numData, int *data, bool scalar)
{
  int numArgs = *numData;
  char buffer[40];
  for (int i = 0; i < numArgs; i++) {
    sprintf(buffer, "%d ", data[i]);
    Tcl_AppendResult(theInterp, buffer, NULL);
  }

  return 0;
}

extern "C" int
OPS_GetDoubleInput(int *numData, double *data)
{
  int size = *numData;
  for (int i = 0; i < size; i++) {
    if ((currentArg >= maxArg) ||
        (Tcl_GetDouble(theInterp, currentArgv[currentArg], &data[i]) !=
         TCL_OK)) {
      // opserr << "OPS_GetDoubleInput -- error reading " << currentArg <<
      // endln;
      return -1;
    } else
      currentArg++;
  }

  return 0;
}

extern "C" int
OPS_SetDoubleOutput(int *numData, double *data, bool scalar)
{
  int numArgs = *numData;
  char buffer[40];
  for (int i = 0; i < numArgs; i++) {
    sprintf(buffer, "%35.20f ", data[i]);
    Tcl_AppendResult(theInterp, buffer, NULL);
  }

  return 0;
}

extern "C" const char *
OPS_GetString(void)
{
  const char *res = 0;
  if (currentArg >= maxArg) {
    // opserr << "OPS_GetStringInput -- error reading " << currentArg << endln;
    return res;
  }
  res = currentArgv[currentArg];

  currentArg++;

  return res;
}

extern "C" const char *
OPS_GetStringFromAll(char *buffer, int len)
{
  return OPS_GetString(); // Everything's a string in Tcl
}

extern "C" int
OPS_SetString(const char *str)
{
  Tcl_SetResult(theInterp, (char *)str, TCL_VOLATILE);
  return 0;
}

int
OPS_GetStringCopy(char **arrayData)
{
  if (currentArg >= maxArg) {
    opserr << "OPS_GetStringInput -- error reading " << currentArg << endln;
    return -1;
  }
  char *newData = new char[strlen(currentArgv[currentArg]) + 1];
  strcpy(newData, currentArgv[currentArg]);
  *arrayData = newData;
  currentArg++;

  return 0;
}

extern "C" matObj *
OPS_GetMaterial(int *matTag, int *matType)
{
  if (*matType == OPS_UNIAXIAL_MATERIAL_TYPE) {
    UniaxialMaterial *theUniaxialMaterial = OPS_getUniaxialMaterial(*matTag);

    if (theUniaxialMaterial != 0) {

      UniaxialMaterial *theCopy = theUniaxialMaterial->getCopy();
      //  uniaxialMaterialObjectCount++;
      // theUniaxialMaterials[uniaxialMaterialObjectCount] = theCopy;

      matObject *theMatObject = new matObject;
      theMatObject->tag = *matTag;
      theMatObject->nParam = 1;
      theMatObject->nState = 0;

      theMatObject->theParam = new double[1];
      //  theMatObject->theParam[0] = uniaxialMaterialObjectCount;
      theMatObject->theParam[0] = 1; // code for uniaxial material

      theMatObject->tState = 0;
      theMatObject->cState = 0;
      theMatObject->matFunctPtr = OPS_InvokeMaterialObject;

      theMatObject->matObjectPtr = theCopy;

      return theMatObject;
    }

    fprintf(stderr, "getMaterial - no uniaxial material exists with tag %d\n",
            *matTag);
    return 0;

  } else if (*matType == OPS_SECTION_TYPE) {
    fprintf(stderr, "getMaterial - not yet implemented for Section\n");
    return 0;
  } else {

    //    NDMaterial *theNDMaterial = theModelBuilder->getNDMaterial(*matTag);

    //    if (theNDMaterial != 0)
    //      theNDMaterial = theNDMaterial->getCopy(matType);
    //    else {
    //      fprintf(stderr,"getMaterial - no nd material exists with tag %d\n",
    //      *matTag); return 0;
    //    }

    //    if (theNDMaterial == 0) {
    //      fprintf(stderr,"getMaterial - material with tag %d cannot deal with
    //      %d\n", *matTag, matType); return 0;
    //    }

    fprintf(stderr, "getMaterial - not yet implemented for nDMaterial\n");
    return 0;
  }

  fprintf(stderr, "getMaterial - unknown material type\n");
  return 0;
}

/*
extern "C"
void OPS_GetMaterialPtr(int *matTag, matObj *theRes)
{
  UniaxialMaterial *theUniaxialMaterial =
theModelBuilder->getUniaxialMaterial(*matTag);

  if (theUniaxialMaterial != 0) {

    UniaxialMaterial *theCopy = theUniaxialMaterial->getCopy();
    if (theCopy  == 0) {
      fprintf(stderr,"OPS_GetMaterialPtr() failed - no material of type %d \n",
*matTag); theRes = 0; return;
    }

    uniaxialMaterialObjectCount++;
    theUniaxialMaterials[uniaxialMaterialObjectCount] = theCopy;

    matObject *theMatObject = new matObject;
    theMatObject->tag = *matTag;
    theMatObject->nParam = 1;
    theMatObject->nState = 0;

    theMatObject->theParam = new double[1];
    theMatObject->theParam[0] = uniaxialMaterialObjectCount;

    theMatObject->tState = 0;
    theMatObject->cState = 0;
    theMatObject->matFunctPtr = OPS_UniaxialMaterialFunction;

    theRes = theMatObject;
  }

  theRes = 0;
}
*/

// 
// END INTERPRETER STUFF
//

extern "C" eleObj *
OPS_GetElement(int *eleTag)
{
  return 0;
}

extern "C" eleObj *
OPS_GetElementType(char *type, int sizeType)
{

  // try existing loaded routines

  ElementFunction *eleFunction = theElementFunctions;
  bool found = false;
  while (eleFunction != NULL && found == false) {
    if (strcmp(type, eleFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      eleObj *theEleObject = new eleObj;
      theEleObject->eleFunctPtr = eleFunction->theFunct;
      return theEleObject;
    } else
      eleFunction = eleFunction->next;
  }

  // ty to load new routine from dynamic library in load path

  eleFunct eleFunctPtr;
  void *libHandle;

  int res = getLibraryFunction(type, type, &libHandle, (void **)&eleFunctPtr);

  if (res == 0) {

    // add the routine to the list of possible elements

    char *funcName = new char[strlen(type) + 1];
    strcpy(funcName, type);
    eleFunction = new ElementFunction;
    eleFunction->theFunct = eleFunctPtr;
    eleFunction->funcName = funcName;
    eleFunction->next = theElementFunctions;
    theElementFunctions = eleFunction;

    // create a new eleObject, set the function ptr &  return it

    eleObj *theEleObject = new eleObj;
    // eleObj *theEleObject = (eleObj *)malloc(sizeof( eleObj));;

    theEleObject->eleFunctPtr = eleFunction->theFunct;

    return theEleObject;
  }

  return 0;
}

extern "C" matObj *
OPS_GetMaterialType(char *type, int sizeType)
{

  // try existing loaded routines
  MaterialFunction *matFunction = theMaterialFunctions;
  bool found = false;
  while (matFunction != NULL && found == false) {
    if (strcmp(type, matFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      matObj *theMatObject = new matObj;
      theMatObject->matFunctPtr = matFunction->theFunct;
      /* opserr << "matObj *OPS_GetMaterialType() - FOUND " << endln;  */
      return theMatObject;
    } else
      matFunction = matFunction->next;
  }

  // ty to load new routine from dynamic library in load path
  matFunct matFunctPtr;
  void *libHandle;

  int res = getLibraryFunction(type, type, &libHandle, (void **)&matFunctPtr);

  if (res == 0) {

    // add the routine to the list of possible elements

    char *funcName = new char[strlen(type) + 1];
    strcpy(funcName, type);
    matFunction = new MaterialFunction;
    matFunction->theFunct = matFunctPtr;
    matFunction->funcName = funcName;
    matFunction->next = theMaterialFunctions;
    theMaterialFunctions = matFunction;

    // create a new eleObject, set the function ptr &  return it

    matObj *theMatObject = new matObj;
    // eleObj *theEleObject = (eleObj *)malloc(sizeof( eleObj));;

    theMatObject->matFunctPtr = matFunction->theFunct;

    //    fprintf(stderr,"getMaterial Address %p\n",theMatObject);

    return theMatObject;
  }

  return 0;
}

extern "C" limCrvObj *
OPS_GetLimitCurveType(char *type, int sizeType)
{

  // try existing loaded routines
  LimitCurveFunction *limCrvFunction = theLimitCurveFunctions;
  bool found = false;
  while (limCrvFunction != NULL && found == false) {
    if (strcmp(type, limCrvFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      limCrvObj *theLimCrvObject = new limCrvObj;
      theLimCrvObject->limCrvFunctPtr = limCrvFunction->theFunct;
      /* opserr << "limCrvObj *OPS_GetLimitCurveType() - FOUND " << endln;  */
      return theLimCrvObject;
    } else
      limCrvFunction = limCrvFunction->next;
  }

  // try to load new routine from dynamic library in load path
  limCrvFunct limCrvFunctPtr;
  void *libHandle;
  int res =
      getLibraryFunction(type, type, &libHandle, (void **)&limCrvFunctPtr);

  if (res == 0) {
    // add the routine to the list of possible elements
    char *funcName = new char[strlen(type) + 1];
    strcpy(funcName, type);
    limCrvFunction = new LimitCurveFunction;
    limCrvFunction->theFunct = limCrvFunctPtr;
    limCrvFunction->funcName = funcName;
    limCrvFunction->next = theLimitCurveFunctions;
    theLimitCurveFunctions = limCrvFunction;

    // create a new eleObject, set the function ptr &  return it
    limCrvObj *theLimCrvObject = new limCrvObj;
    theLimCrvObject->limCrvFunctPtr = limCrvFunction->theFunct;
    return theLimCrvObject;
  }

  return 0;
}

extern "C" int
OPS_AllocateLimitCurve(limCrvObject *theLimCrv)
{

  /*fprintf(stderr,"allocateLimitCurve Address %p\n",theLimCrv);*/

  if (theLimCrv->nParam > 0)
    theLimCrv->theParam = new double[theLimCrv->nParam];

  int nState = theLimCrv->nState;

  if (nState > 0) {
    theLimCrv->cState = new double[nState];
    theLimCrv->tState = new double[nState];
    for (int i = 0; i < nState; i++) {
      theLimCrv->cState[i] = 0;
      theLimCrv->tState[i] = 0;
    }
  } else {
    theLimCrv->cState = 0;
    theLimCrv->tState = 0;
  }

  return 0;
}

extern "C" int
OPS_AllocateMaterial(matObject *theMat)
{

  /*fprintf(stderr,"allocateMaterial Address %p\n",theMat);*/

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
    /*  opserr << "AllocateElement - matTag " << matTags[i] << "\n"; */

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

UniaxialMaterial *
Tcl_addWrapperUniaxialMaterial(matObj *theMat, ClientData clientData,
                               Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theInterp = interp;

  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

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
    return 0;
  }

  WrapperUniaxialMaterial *theMaterial =
      new WrapperUniaxialMaterial(argv[1], theMat);

  return theMaterial;
}

NDMaterial *
Tcl_addWrapperNDMaterial(matObj *theMat, ClientData clientData,
                         Tcl_Interp *interp, int argc, TCL_Char **argv,
                         TclBasicBuilder *theModelbuilder)
{
  theInterp = interp;

  // theModelBuilder = builder;
  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

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

  if (result != 0 ||
      (matType != OPS_PLANESTRESS_TYPE && matType != OPS_PLANESTRAIN_TYPE &&
       matType != OPS_THREEDIMENSIONAL_TYPE)) {
    opserr << "Tcl_addWrapperNDMaterial - failed in element function " << result
           << endln;
    return 0;
  }

  WrapperNDMaterial *theMaterial =
      new WrapperNDMaterial(argv[1], theMat, theMat->matType);

  return theMaterial;
}

LimitCurve *
Tcl_addWrapperLimitCurve(limCrvObj *theLimCrv, ClientData clientData,
                         Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theInterp = interp;

  //  theModelBuilder = builder;
  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

  // get the current load factor
  static modelState theModelState;
  if (theDomain != 0) {
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;
    theModelState.time = time;
    theModelState.dt = dt;
  }

  // invoke the limit curve function with isw = 0
  int isw = ISW_INIT;
  int result;
  theLimCrv->limCrvFunctPtr(theLimCrv, &theModelState, 0, 0, 0, &isw, &result);

  if (result != 0) {
    opserr << "Tcl_addWrapperLimitCurve - failed in limit curve function "
           << result << endln;
    return 0;
  }

  WrapperLimitCurve *theLimitCurve = new WrapperLimitCurve(argv[1], theLimCrv);

  return theLimitCurve;
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

int
G3_raise(G3_Runtime *rt, const char *msg, ...){
  va_list ap;

  va_start(ap, msg);
  int n = vsnprintf(NULL, 0, msg, ap);
  va_end(ap);

  if (n < 0)
    return -1;

  size_t size = (size_t)n + 1 + 8;
  char *new_str = (char*)malloc(size);
  if (new_str == NULL)
    return -1;

  strcpy(new_str, "error {");
  va_start(ap, msg);
  n = vsnprintf(new_str+7, size, msg, ap);
  va_end(ap);
  strcpy(new_str+7+n, "}\n");


  Tcl_Interp *tcl_interp = G3_getInterpreter(rt);
  Tcl_Eval(tcl_interp, new_str);
  Tcl_Obj *infoObj = Tcl_GetVar2Ex(tcl_interp, "errorInfo", NULL, TCL_GLOBAL_ONLY);
  const char * error_str = Tcl_GetString(infoObj);
  opserr << error_str;

  /*
  Tcl_Obj *top_interpInfoName ;
  Tcl_Obj *top_interpInfo ;
    top_interpInfoName = Tcl_NewStringObj("errorInfo", -1) ;
    Tcl_IncrRefCount(top_interpInfoName) ;
    top_interpInfo =  Tcl_ObjGetVar2(tcl_interp,
                                     top_interpInfoName,
                                     NULL,
                                     TCL_LEAVE_ERR_MSG) ;
    Tcl_IncrRefCount(top_interpInfo) ;
    const char *error_str = Tcl_GetString(top_interpInfo);
    opserr << "ERROR -- " << msg << "\n\n" << error_str;
    Tcl_DecrRefCount(top_interpInfoName) ;
    Tcl_DecrRefCount(top_interpInfo);
    */
    // throw 20;
    return 0;
}

G3_Runtime *
G3_getRuntime(Tcl_Interp *interp)
{
  G3_Runtime *rt = (G3_Runtime*)Tcl_GetAssocData(interp, "G3_Runtime", NULL);
  if (!rt)
    opserr << G3_WARN_PROMPT << " No runtime\n";;
  return rt;
}

Tcl_Interp *
G3_getInterpreter(G3_Runtime* rt)
{return rt->m_interp;}

TclBuilder *
G3_getModelBuilder(G3_Runtime *rt)
{return rt->m_builder;}

int
G3_setModelBuilder(G3_Runtime *rt, TclBuilder* builder)
{
  rt->m_builder = builder;
  return 1;
}

TclSafeBuilder *
G3_getSafeBuilder(G3_Runtime *rt)
{
  return (TclSafeBuilder*)G3_getModelBuilder(rt);
  /*
  Tcl_Interp *interp = G3_getInterpreter(rt);
  TclSafeBuilder *theTclBuilder =
      (TclSafeBuilder *)Tcl_GetAssocData(interp, "OPS::theTclSafeBuilder", NULL);
  return theTclBuilder;
  */
}

int
G3_setDomain(G3_Runtime *rt, Domain* domain){
  int exists = rt->m_domain ? 1 : 0;
  Domain *old = rt->m_domain;
  /*
  if (old && old != domain)
    throw 20;
  printf("SETTING: %p->%p\n", rt, rt->m_domain);
  */
  rt->m_domain = domain;
  // opserr << "Domain set from '" << (long int)old << "' to '" << (long int)domain << "'\n";
  return exists;
}


Domain *
G3_getDomain(G3_Runtime *rt)
{
  // Tcl_Interp *interp = G3_getInterpreter(rt);
  return rt->m_domain;
}

int G3_addTimeSeries(G3_Runtime *rt, TimeSeries *series)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  TclSafeBuilder *builder = G3_getSafeBuilder(rt);
      // (TclSafeBuilder *)Tcl_GetAssocData(interp, "OPS::theTclSafeBuilder", NULL);
  return builder->addTimeSeries(series);
}

int G3_removeTimeSeries(G3_Runtime *rt, int tag) {
  // Tcl_Interp *interp = G3_getInterpreter(rt);
  TclSafeBuilder *builder = G3_getSafeBuilder(rt);
  if (builder)
    // TODO
    return true;
  else
    return false;
}


TimeSeries *G3_getTimeSeries(G3_Runtime *rt, int tag)
{
  // Tcl_Interp *interp = G3_getInterpreter(rt);

  TimeSeries *series;
  // TclSafeBuilder *builder =
  TclSafeBuilder *builder = G3_getSafeBuilder(rt);
      // (TclSafeBuilder *)Tcl_GetAssocData(interp, "OPS::theTclSafeBuilder", NULL);
  if (builder) {
     series = builder->getTimeSeries(std::to_string(tag));
  } else {
    // opserr << "WARNING: Unable to find safe model builder\n";
    series = nullptr;
  }

  return series;
}

/*
void G3_clearAllTimeSeries(void) {
  theTimeSeriesObjects.clearAll();
}
*/   


extern "C" int
OPS_InvokeMaterialDirectly2(matObject *theMat, modelState *model,
                            double *strain, double *stress, double *tang,
                            int *isw)
{
  int error = 0;
  //  fprintf(stderr,"invokeMaterialDirectly Address %d %d\n",theMat,
  //  sizeof(int));
  if (theMat != 0)
    theMat->matFunctPtr(theMat, model, strain, tang, stress, isw, &error);
  else
    error = -1;

  return error;
}

#if !defined(OPS_USE_RUNTIME)
UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag)
{
  return OPS_getUniaxialMaterial(matTag);
}
#endif

CrdTransf *
G3_getCrdTransf(G3_Runtime *rt, G3_Tag tag)
{
  TclSafeBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    return nullptr;
  }
  return builder->getCrdTransf(tag);
}

UniaxialMaterial *
G3_getUniaxialMaterialInstance(G3_Runtime *rt, int tag)
{
  TclSafeBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    // TODO
    return OPS_getUniaxialMaterial(tag);
  }

  UniaxialMaterial *mat = builder->getUniaxialMaterial(tag);
  // TODO
  return mat ? mat : OPS_getUniaxialMaterial(tag);
}

int G3_addUniaxialMaterial(G3_Runtime *rt, UniaxialMaterial *mat) {
  TclSafeBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    opserr << "WARNING Failed to find safe model builder\n";
    return 0;
  }
  int stat = builder->addUniaxialMaterial(mat);
  return stat ? TCL_OK : TCL_ERROR;
}


NDMaterial *
OPS_GetNDMaterial(int matTag)
{
  return OPS_getNDMaterial(matTag);
}

SectionForceDeformation *
OPS_GetSectionForceDeformation(int secTag)
{
  return OPS_getSectionForceDeformation(secTag);
}

FrictionModel *
OPS_GetFrictionModel(int frnTag) {return OPS_getFrictionModel(frnTag);}

int
OPS_GetNDF() {return theModelBuilder->getNDF();}

bool
G3_modelIsBuilt(G3_Runtime* rt) {return rt->model_is_built;}

int
G3_getNDM(G3_Runtime *rt)
{
  TclBuilder *builder;
  if (builder = G3_getModelBuilder(rt))
    return builder->getNDM();
  else
    return -1;
}

int
G3_getNDF(G3_Runtime *rt)
{
  TclBuilder *builder;
  if (builder = G3_getModelBuilder(rt))
    return builder->getNDF();
  else
    return -1;
}

int
OPS_GetNDM(void) {return theModelBuilder->getNDM();}

FE_Datastore *
OPS_GetFEDatastore() {return theDatabase;}

const char *
OPS_GetInterpPWD() {return getInterpPWD(theInterp);}

#if !defined(OPS_USE_RUNTIME)
  Domain *
  OPS_GetDomain(void) {return theDomain;}

  AnalysisModel **
  OPS_GetAnalysisModel(void){return &theAnalysisModel;}

  CrdTransf *
  OPS_GetCrdTransf(int crdTag) {return OPS_getCrdTransf(crdTag);}

#endif

void
TCL_OPS_setModelBuilder(TclBasicBuilder *theNewBuilder) {theModelBuilder = theNewBuilder;}

LimitCurve *
OPS_GetLimitCurve(int LimCrvTag)
{
  return OPS_getLimitCurve(LimCrvTag);
}

EquiSolnAlgo **
OPS_GetAlgorithm(void)
{
  return &theAlgorithm;
}

ConstraintHandler **
OPS_GetHandler(void)
{
  return &theHandler;
}

DOF_Numberer **
OPS_GetNumberer(void)
{
  return &theNumberer;
}

LinearSOE **
OPS_GetSOE(void)
{
  return &theSOE;
}

LinearSOE **
G3_getLinearSoePtr(G3_Runtime* rt) {
  LinearSOE** soe =  &rt->m_sys_of_eqn;
  // opserr << "DEBUG G3_getLinearSoe(" << (void*)rt << ") -> " << (void*)(*soe) << "\n";
  return soe;
}

int G3_setLinearSoe(G3_Runtime* rt, LinearSOE* soe)
{
  // opserr << "DEBUG G3_setLinearSoe(" << (void*)rt << (void*)soe << ")\n";
  rt->m_sys_of_eqn = soe;
  // if the analysis exists - we want to change the SOE
  if (soe != nullptr) {
    StaticAnalysis* static_analysis = G3_getStaticAnalysis(rt);
    if (static_analysis != 0)
      static_analysis->setLinearSOE(*soe);
    
    DirectIntegrationAnalysis *direct_trans_analysis = G3_getTransientAnalysis(rt);
    if (direct_trans_analysis != 0)
      direct_trans_analysis->setLinearSOE(*soe);

#ifdef _PARALLEL_PROCESSING
      if (static_analysis != 0 || direct_trans_analysis != 0) {
        SubdomainIter &theSubdomains = theDomain.getSubdomains();
        Subdomain *theSub;
        while ((theSub = theSubdomains()) != 0) {
          theSub->setAnalysisLinearSOE(*theSOE);
        }
      }
#endif
  }
  return 0;
}

LinearSOE *
G3_getDefaultLinearSoe(G3_Runtime* rt, int flags) {
  // `flags` is unused right now but could be useful
  // for ensuring properties about the SOE, like
  // forcing fullGen.
  LinearSOE* theSOE = *G3_getLinearSoePtr(rt);

  opsdbg << "DEBUG G3_getDefaultLinearSoe(" << (long int)rt << ", " << flags << ")-> " << (long int)theSOE << "\n";

  if (theSOE == NULL) {
    opswrn << "no LinearSOE specified, default ProfileSPDLinSOE will be used\n";
    ProfileSPDLinSolver *theSolver;
    theSolver = new ProfileSPDLinDirectSolver();
#ifdef _PARALLEL_PROCESSING
    theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
    theSOE = new ProfileSPDLinSOE(*theSolver);
#endif
    G3_setLinearSoe(rt, theSOE);
  }
  return theSOE;
}



EigenSOE **
OPS_GetEigenSOE(void)
{
  return &theEigenSOE;
}

StaticAnalysis **
OPS_GetStaticAnalysis(void)
{
  return &theStaticAnalysis;
}

int
G3_setAnalysisModel(G3_Runtime *rt, AnalysisModel *the_analysis)
{
  rt->m_analysis_model = the_analysis;
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theAnalysisModel", NULL, (ClientData)the_analysis);
  return 1;
}

AnalysisModel *
G3_getAnalysisModel(G3_Runtime *rt){return rt->m_analysis_model;}

AnalysisModel **
G3_getAnalysisModelPtr(G3_Runtime *rt){return rt->m_analysis_model_ptr;}


StaticAnalysis *
G3_getStaticAnalysis(G3_Runtime *rt)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  StaticAnalysis *analysis =
      (StaticAnalysis *)Tcl_GetAssocData(interp, "OPS::theStaticAnalysis", NULL);
  return analysis;
}


int
G3_setStaticAnalysis(G3_Runtime *rt, StaticAnalysis *the_analysis)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theStaticAnalysis", NULL, (ClientData)the_analysis);
  return 1;
}

int
G3_delStaticAnalysis(G3_Runtime *rt)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  StaticAnalysis* ana;
  // TODO
  if (ana=G3_getStaticAnalysis(rt))
    ;// delete ana;
  Tcl_SetAssocData(interp, "OPS::theStaticAnalysis", NULL, (ClientData)NULL);
  return 1;
}

StaticIntegrator *
G3_getStaticIntegrator(G3_Runtime *rt)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  StaticIntegrator *tsi =
      (StaticIntegrator *)Tcl_GetAssocData(interp, "OPS::theStaticIntegrator", NULL);

  return tsi;
}

int
G3_setStaticIntegrator(G3_Runtime *rt, StaticIntegrator *the_analysis)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theStaticIntegrator", NULL, (ClientData)the_analysis);
  return 1;
}


DirectIntegrationAnalysis **
OPS_GetTransientAnalysis(void)
{
  return &theTransientAnalysis;
}

DirectIntegrationAnalysis *
G3_getTransientAnalysis(G3_Runtime *rt)
{
  // TODO
  /*
  Tcl_Interp *interp = G3_getInterpreter(rt);
  DirectIntegrationAnalysis *analysis =
      (DirectIntegrationAnalysis *)Tcl_GetAssocData(interp, "OPS::theTransientAnalysis", NULL);
  */
  return theTransientAnalysis;

}

int
G3_setTransientAnalysis(G3_Runtime *rt, DirectIntegrationAnalysis *the_analysis)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theTransientAnalysis", NULL, (ClientData)the_analysis);
  theTransientAnalysis = the_analysis;
  return 1;
}


VariableTimeStepDirectIntegrationAnalysis **
OPS_GetVariableTimeStepTransientAnalysis(void)
{
  return &theVariableTimeStepTransientAnalysis;
}

int *
OPS_GetNumEigen(void)
{
  return &numEigen;
}

StaticIntegrator **
OPS_GetStaticIntegrator(void)
{
  return &theStaticIntegrator;
}

TransientIntegrator **
OPS_GetTransientIntegrator(void)
{
  return &theTransientIntegrator;
}

ConvergenceTest **
OPS_GetTest(void)
{
  return &theTest;
}

bool *
OPS_builtModel(void)
{
  return &builtModel;
}

int
OPS_numIter()
{
  return 0;
}
