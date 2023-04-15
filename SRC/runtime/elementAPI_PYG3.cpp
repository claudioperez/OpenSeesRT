/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <assert.h>
#include <elementAPI.h>
#include <stdlib.h>
#include <packages.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <g3_api.h>
#include <G3_Runtime.h>
#include <G3_Logging.h>
#include <runtime/BasicModelBuilder.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <FrictionModel.h>

#include <DirectIntegrationAnalysis.h>
#include <StaticAnalysis.h>

// extern AnalysisModel             *theAnalysisModel;
// extern ConstraintHandler         *theHandler;
// extern DOF_Numberer              *theGlobalNumberer;
// extern ConvergenceTest           *theTest;
extern EquiSolnAlgo              *theAlgorithm;
extern LinearSOE                 *theSOE;
extern EigenSOE                  *theEigenSOE;
extern StaticAnalysis            *theStaticAnalysis;
extern StaticIntegrator          *theStaticIntegrator;
extern DirectIntegrationAnalysis *theTransientAnalysis;
extern TransientIntegrator       *theTransientIntegrator;
extern VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;
extern bool builtModel;
extern FE_Datastore *theDatabase;

static Tcl_Interp *theInterp       = nullptr;
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
OPS_Error(const char *errorMessage, int length)
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
               TCL_Char ** const argv, Domain *domain, TclBuilder *builder)
{
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;
  return 0;
}


extern "C" int
OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                        int mArg, TCL_Char ** const argv, Domain *domain)
{
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;
  return 0;
}


#if 0
extern "C" int
G3_GetIntInput(Tcl_Interp* interp, const char**const argv, int* argc, int *numData, int *data)
{
  int size = *numData;

  for (int i = 0; i < size; i++) {
    if ((currentArg >= maxArg) || (Tcl_GetInt(theInterp, argv[argc], &data[i]) != TCL_OK)) {
      // opserr << "OPS_GetIntInput -- error reading " << currentArg << endln;
      return -1;
    } else
      currentArg++;
  }

  return 0;
}
#endif

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

G3_Runtime *
G3_getRuntime(Tcl_Interp *interp)
{
  G3_Runtime *rt = (G3_Runtime*)Tcl_GetAssocData(interp, "G3_Runtime", nullptr);
  if (!rt)
    opserr << G3_WARN_PROMPT << " No runtime\n";;
  return rt;
}

Tcl_Interp *
G3_getInterpreter(G3_Runtime* rt) {return rt->m_interp;}

TclBuilder *
G3_getModelBuilder(G3_Runtime *rt) {return rt->m_builder;}

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
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
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


CrdTransf *
G3_getCrdTransf(G3_Runtime *rt, G3_Tag tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    return nullptr;
  }
  return builder->getCrdTransf(tag);
}

#if 0
SectionForceDeformation *
OPS_GetSectionForceDeformation(int secTag)
{
  return OPS_getSectionForceDeformation(secTag);
}
#endif

SectionForceDeformation*
G3_getSectionForceDeformation(G3_Runtime* rt, int tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  assert(builder);
  return builder->getSection(tag);
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
G3_GetNDMaterial(G3_Runtime* rt, int matTag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  return builder->getNDMaterial(matTag);
}

NDMaterial *
OPS_GetNDMaterial(int matTag)
{
  return OPS_getNDMaterial(matTag);
}


FrictionModel *
OPS_GetFrictionModel(int frnTag) {return OPS_getFrictionModel(frnTag);}

int
OPS_GetNDF() {return theModelBuilder->getNDF();}

int
G3_getNDM(G3_Runtime *rt)
{
  BasicModelBuilder *builder;
  if ((builder = G3_getSafeBuilder(rt)))
    return builder->getNDM();
  else
    return -1;
}

int
OPS_GetNDM(void) {return theModelBuilder->getNDM();}

bool *
OPS_builtModel(void) {return &builtModel;}

LinearSOE **
G3_getLinearSoePtr(G3_Runtime* rt) {
  LinearSOE** soe =  &rt->m_sys_of_eqn;
  return soe;
}


AnalysisModel **
G3_getAnalysisModelPtr(G3_Runtime *rt){return rt->m_analysis_model_ptr;}

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

FE_Datastore *
OPS_GetFEDatastore() {return theDatabase;}

const char *
OPS_GetInterpPWD() {return getInterpPWD(theInterp);}

EquiSolnAlgo **
OPS_GetAlgorithm(void) {return &theAlgorithm;}

LinearSOE **
OPS_GetSOE(void) {return &theSOE;}

EigenSOE **
OPS_GetEigenSOE(void) {return &theEigenSOE;}

StaticAnalysis **
OPS_GetStaticAnalysis(void) {return &theStaticAnalysis;}

#if 0 && !defined(OPS_USE_RUNTIME)
UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag)
{
  return OPS_getUniaxialMaterial(matTag);
}

Domain *
OPS_GetDomain(void) {return theDomain;}

AnalysisModel **
OPS_GetAnalysisModel(void){return &theAnalysisModel;}

CrdTransf *
OPS_GetCrdTransf(int crdTag) {return OPS_getCrdTransf(crdTag);}

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

