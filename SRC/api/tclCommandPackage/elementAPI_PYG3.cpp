/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

/*
 * Written: cmp
 */

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

extern AnalysisModel            *theAnalysisModel;
extern EquiSolnAlgo             *theAlgorithm;
extern ConstraintHandler        *theHandler;
extern DOF_Numberer             *theGlobalNumberer;
extern LinearSOE                *theSOE;
extern EigenSOE                 *theEigenSOE;
extern StaticAnalysis           *theStaticAnalysis;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;
extern StaticIntegrator *theStaticIntegrator;
extern TransientIntegrator *theTransientIntegrator;
extern ConvergenceTest *theTest;
extern bool builtModel;
extern FE_Datastore *theDatabase;

static ElementFunction  *theElementFunctions = NULL;
static MaterialFunction *theMaterialFunctions = NULL;

static Tcl_Interp *theInterp = 0;
static Domain     *theDomain = 0;
static TclBuilder *theModelBuilder = 0;

static TCL_Char **currentArgv = 0;
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
//  G3_setDomain(rt, domain);
//  G3_setModelBuilder(rt, builder);
// theInterp = interp;
// theDomain = domain;
// theModelBuilder = builder;
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
//  G3_setDomain(rt, domain);
//  theInterp = interp;
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
/*
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
*/

/*
extern "C"
void OPS_GetMaterialPtr(int *matTag, matObj *theRes)
{
  UniaxialMaterial *theUniaxialMaterial = theModelBuilder->getUniaxialMaterial(*matTag);

  if (theUniaxialMaterial != 0) {

    UniaxialMaterial *theCopy = theUniaxialMaterial->getCopy();
    if (theCopy  == 0) {
      fprintf(stderr,"OPS_GetMaterialPtr() failed - no material of type %d \n", *matTag); 
      theRes = 0; 
      return;
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
OPS_GetElement(int *eleTag) {return 0;}

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
  /*
  Tcl_Interp *interp = G3_getInterpreter(rt);
  BasicModelBuilder *theTclBuilder =
      (BasicModelBuilder *)Tcl_GetAssocData(interp, "OPS::theBasicModelBuilder", NULL);
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
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
      // (BasicModelBuilder *)Tcl_GetAssocData(interp, "OPS::theBasicModelBuilder", NULL);
  return builder->addTimeSeries(series);
}

/*
int G3_removeTimeSeries(G3_Runtime *rt, int tag) {
  // Tcl_Interp *interp = G3_getInterpreter(rt);
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
  if (builder)
    // TODO
    return true;
  else
    return false;
}
*/

TimeSeries *G3_getTimeSeries(G3_Runtime *rt, int tag)
{
  // Tcl_Interp *interp = G3_getInterpreter(rt);

  TimeSeries *series;
  // BasicModelBuilder *builder =
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
      // (BasicModelBuilder *)Tcl_GetAssocData(interp, "OPS::theBasicModelBuilder", NULL);
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

void
TCL_OPS_setModelBuilder(TclBasicBuilder *theNewBuilder) {theModelBuilder = theNewBuilder;}

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


/*
DirectIntegrationAnalysis **
OPS_GetTransientAnalysis(void) {return &theTransientAnalysis;}

DirectIntegrationAnalysis *
G3_getTransientAnalysis(G3_Runtime *rt)
{
  // TODO
  
  // Tcl_Interp *interp = G3_getInterpreter(rt);
  // DirectIntegrationAnalysis *analysis =
  //     (DirectIntegrationAnalysis *)Tcl_GetAssocData(interp, "OPS::theTransientAnalysis", NULL);
  // 
  return theTransientAnalysis;

}
*/

StaticIntegrator **
OPS_GetStaticIntegrator(void) {return &theStaticIntegrator;}

TransientIntegrator **
OPS_GetTransientIntegrator(void) {return &theTransientIntegrator;}

/*
ConvergenceTest **
OPS_GetTest(void) {return &theTest;}

ConstraintHandler **
OPS_GetHandler(void) {return &theHandler;}

DOF_Numberer **
OPS_GetNumberer(void) {return &theGlobalNumberer;}

int G3_setLinearSoe(G3_Runtime* rt, LinearSOE* soe)
{
  rt->m_sys_of_eqn = soe;
  // if the analysis exists - we want to change the SOE
  if (soe != nullptr) {
    StaticAnalysis* static_analysis = G3_getStaticAnalysis(rt);
    if (static_analysis != 0)
      static_analysis->setLinearSOE(*soe);
    
    DirectIntegrationAnalysis *direct_trans_analysis = theTransientAnalysis; // G3_getTransientAnalysis(rt);
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

bool
G3_modelIsBuilt(G3_Runtime* rt) {return rt->model_is_built;}


*/


bool *
OPS_builtModel(void) {return &builtModel;}


