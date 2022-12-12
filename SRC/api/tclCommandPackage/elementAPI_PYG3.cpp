/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
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

#include <runtime/BasicModelBuilder.h>
#include <runtime/BasicModelBuilder.h>
#include <WrapperElement.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <FrictionModel.h>

#include <DirectIntegrationAnalysis.h>
#include <StaticAnalysis.h>

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

extern AnalysisModel             *theAnalysisModel;
extern EquiSolnAlgo              *theAlgorithm;
extern ConstraintHandler         *theHandler;
extern DOF_Numberer              *theGlobalNumberer;
extern LinearSOE                 *theSOE;
extern EigenSOE                  *theEigenSOE;
extern StaticAnalysis            *theStaticAnalysis;
extern StaticIntegrator          *theStaticIntegrator;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern TransientIntegrator       *theTransientIntegrator;
extern VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;
extern ConvergenceTest *theTest;
extern bool builtModel;
extern FE_Datastore *theDatabase;

static ElementFunction  *theElementFunctions  = nullptr;
static MaterialFunction *theMaterialFunctions = nullptr;

static Tcl_Interp *theInterp       = nullptr;
static Domain     *theDomain       = nullptr;
static TclBuilder *theModelBuilder = nullptr;

static TCL_Char **currentArgv      = nullptr;
static int currentArg = 0;
static int maxArg = 0;

modelState theModelState;

extern const char *getInterpPWD(Tcl_Interp *interp);

// static int uniaxialMaterialObjectCount = 0;


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
    Tcl_AppendResult(theInterp, buffer, nullptr);
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
    Tcl_AppendResult(theInterp, buffer, nullptr);
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
  return OPS_GetString();
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

// 
// END INTERPRETER STUFF
//

extern "C" eleObj *
OPS_GetElement(int *eleTag) {return 0;}

extern "C" eleObj *
OPS_GetElementType(char *type, int sizeType)
{

  // try existing loaded routines

  ElementFunction *eleFunction = theElementFunctions;
  bool found = false;
  while (eleFunction != nullptr && found == false) {
    if (strcmp(type, eleFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      eleObj *theEleObject = new eleObj;
      theEleObject->eleFunctPtr = eleFunction->theFunct;
      return theEleObject;
    } else
      eleFunction = eleFunction->next;
  }

  // try to load new routine from dynamic library in load path

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

    theEleObject->eleFunctPtr = eleFunction->theFunct;

    return theEleObject;
  }

  return 0;
}

// 
extern "C" matObj *
OPS_GetMaterialType(char *type, int sizeType)
{

  // try existing loaded routines
  MaterialFunction *matFunction = theMaterialFunctions;
  bool found = false;
  while (matFunction != nullptr && found == false) {
    if (strcmp(type, matFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      matObj *theMatObject = new matObj;
      theMatObject->matFunctPtr = matFunction->theFunct;
      /* opserr << "matObj *OPS_GetMaterialType() - FOUND " << endln;  */
      return theMatObject;
    } else
      matFunction = matFunction->next;
  }

  // try to load new routine from dynamic library in load path
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

    // create a new matObject, set the function ptr &  return it

    matObj *theMatObject = new matObj;

    theMatObject->matFunctPtr = matFunction->theFunct;

    return theMatObject;
  }

  return 0;
}

G3_Runtime *
G3_getRuntime(Tcl_Interp *interp)
{
  G3_Runtime *rt = (G3_Runtime*)Tcl_GetAssocData(interp, "G3_Runtime", nullptr);
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
G3_setModelBuilder(G3_Runtime *rt, BasicModelBuilder* builder)
{
  theModelBuilder = builder;
  rt->m_builder = builder;
  return 1;
}

BasicModelBuilder *
G3_getSafeBuilder(G3_Runtime *rt)
{
  return (BasicModelBuilder*)G3_getModelBuilder(rt);
}

int
G3_setDomain(G3_Runtime *rt, Domain* domain){
  int exists = rt->m_domain ? 1 : 0;
  Domain *old = rt->m_domain;
  rt->m_domain = domain;
  return exists;
}


Domain *
G3_getDomain(G3_Runtime *rt)
{
  return rt->m_domain;
}

int G3_addTimeSeries(G3_Runtime *rt, TimeSeries *series)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
      // (BasicModelBuilder *)Tcl_GetAssocData(interp, "OPS::theBasicModelBuilder", NULL);
  return builder->addTimeSeries(series);
}


TimeSeries *G3_getTimeSeries(G3_Runtime *rt, int tag)
{

  TimeSeries *series;
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
  if (builder) {
     series = builder->getTimeSeries(std::to_string(tag));
  } else {
    series = nullptr;
  }

  return series;
}


extern "C" int
OPS_InvokeMaterialDirectly2(matObject *theMat, modelState *model,
                            double *strain, double *stress, double *tang,
                            int *isw)
{
  int error = 0;
  if (theMat != nullptr)
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
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    return nullptr;
  }
  return builder->getCrdTransf(tag);
}

UniaxialMaterial *
G3_getUniaxialMaterialInstance(G3_Runtime *rt, int tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    // TODO
    return OPS_getUniaxialMaterial(tag);
  }

  UniaxialMaterial *mat = builder->getUniaxialMaterial(tag);
  // TODO
  return mat ? mat : OPS_getUniaxialMaterial(tag);
}

int G3_addUniaxialMaterial(G3_Runtime *rt, UniaxialMaterial *mat) {
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    opserr << "WARNING Failed to find safe model builder\n";
    return 0;
  }
  return builder->addUniaxialMaterial(mat);
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

int
G3_getNDM(G3_Runtime *rt)
{
  BasicModelBuilder *builder;
  if (builder = G3_getSafeBuilder(rt))
    return builder->getNDM();
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

EquiSolnAlgo **
OPS_GetAlgorithm(void) {return &theAlgorithm;}


LinearSOE **
OPS_GetSOE(void) {return &theSOE;}

LinearSOE **
G3_getLinearSoePtr(G3_Runtime* rt) {
  LinearSOE** soe =  &rt->m_sys_of_eqn;
  return soe;
}


EigenSOE **
OPS_GetEigenSOE(void) {return &theEigenSOE;}

StaticAnalysis **
OPS_GetStaticAnalysis(void) {return &theStaticAnalysis;}

bool *
OPS_builtModel(void) {return &builtModel;}

#if 0
int
G3_setAnalysisModel(G3_Runtime *rt, AnalysisModel *the_analysis)
{
  rt->m_analysis_model = the_analysis;
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theAnalysisModel", nullptr, (ClientData)the_analysis);
  return 1;
}

AnalysisModel *
G3_getAnalysisModel(G3_Runtime *rt){return rt->m_analysis_model;}
#endif

AnalysisModel **
G3_getAnalysisModelPtr(G3_Runtime *rt){return rt->m_analysis_model_ptr;}


StaticAnalysis *
G3_getStaticAnalysis(G3_Runtime *rt)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  StaticAnalysis *analysis =
      (StaticAnalysis *)Tcl_GetAssocData(interp, "OPS::theStaticAnalysis", nullptr);
  return analysis;
}


int
G3_setStaticAnalysis(G3_Runtime *rt, StaticAnalysis *the_analysis)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);
  Tcl_SetAssocData(interp, "OPS::theStaticAnalysis", nullptr, (ClientData)the_analysis);
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
  Tcl_SetAssocData(interp, "OPS::theStaticAnalysis", nullptr, nullptr);
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


#if 0
StaticIntegrator **
OPS_GetStaticIntegrator(void) {return &theStaticIntegrator;}

TransientIntegrator **
OPS_GetTransientIntegrator(void) {return &theTransientIntegrator;}

DirectIntegrationAnalysis **
OPS_GetTransientAnalysis(void) {return &theTransientAnalysis;}

ConvergenceTest **
OPS_GetTest(void) {return &theTest;}

ConstraintHandler **
OPS_GetHandler(void) {return &theHandler;}

DOF_Numberer **
OPS_GetNumberer(void) {return &theGlobalNumberer;}
#endif



