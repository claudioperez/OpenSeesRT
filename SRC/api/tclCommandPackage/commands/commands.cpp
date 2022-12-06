/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Written: cmp

// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified.

#include <assert.h>
#include <g3_api.h>
#include <G3_Logging.h>

#include <classTags.h>
#include <DOF_Group.h>

extern "C" {
#include <g3_api.h>
}

#include <OPS_Globals.h>
#include <Matrix.h>
#include <set>
#include <vector>
#include <algorithm>

#ifdef _USING_STL_STREAMS
#  include <iomanip>
#else
#  include <DummyStream.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <elementAPI.h>
#include <g3_api.h>

#include <packages.h>
#include <TclPackageClassBroker.h>

#include <ModelBuilder.h>
#include "commands.h"

// domain
#ifdef _PARALLEL_PROCESSING
#  include <PartitionedDomain.h>
#else
#  include <Domain.h>
#endif

#include <Information.h>
#include <Element.h>
#include <Node.h>
#include <ElementIter.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <ParameterIter.h>
#include <SP_Constraint.h>     //Joey UC Davis
#include <SP_ConstraintIter.h> //Joey UC Davis
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>
#include <Pressure_Constraint.h>

// analysis
#include <AnalysisModel.h>
#include <EquiSolnAlgo.h>
#include <Integrator.h>
#include <StaticIntegrator.h>
#include <Newmark.h>

#include "analysis/analysis.h"

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

#include <EigenSOE.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>

#include <LinearSOE.h>


//  recorders
#include <Recorder.h>

#include <ErrorHandler.h>
#include <ConsoleErrorHandler.h>

#include <Response.h>

//
// Global variables
//
Domain theDomain;
// Domain *theGlobalDomainPtr;

ModelBuilder        *theBuilder = nullptr;

EquiSolnAlgo        *theAlgorithm = nullptr;
ConstraintHandler   *theHandler = nullptr;
DOF_Numberer        *theGlobalNumberer = nullptr;
LinearSOE           *theSOE = nullptr;
EigenSOE            *theEigenSOE = nullptr;
StaticIntegrator    *theStaticIntegrator = nullptr;
TransientIntegrator *theTransientIntegrator = nullptr;

AnalysisModel       *theAnalysisModel = nullptr;
StaticAnalysis      *theStaticAnalysis = nullptr;
DirectIntegrationAnalysis *theTransientAnalysis = nullptr;
VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis = nullptr;

ConvergenceTest *theTest = nullptr;
bool builtModel = false;

static char *resDataPtr = nullptr;
static int resDataSize = 0;

#include <FileStream.h>

TclPackageClassBroker theBroker;



// extern int G3_AddTclAnalysisAPI(Tcl_Interp *interp, Domain* domain);
const char *getInterpPWD(Tcl_Interp *interp);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char **argv, Domain *domain);
int printModelGID(ClientData, Tcl_Interp *, int, TCL_Char **);

Tcl_CmdProc TclCommand_record;
Tcl_CmdProc maxOpenFiles;
Tcl_CmdProc setPrecision;
Tcl_CmdProc logFile;
Tcl_CmdProc version;


// TODO: reimplement  int defaultUnits(ClientData, Tcl_Interp *, int, TCL_Char **argv);
// int setParameter(ClientData, Tcl_Interp *, int, TCL_Char **);

// extern
int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

Tcl_CmdProc TclCommand_setLoadConst;
Tcl_CmdProc TclCommand_setCreep;

// domain.cpp
Tcl_CmdProc domainChange;


int
G3_AddTclDomainCommands(Tcl_Interp *interp, Domain* the_domain)
{

  ClientData domain = (ClientData)the_domain;

  Tcl_CreateCommand(interp, "algorithmRecorder", &addAlgoRecorder, domain, nullptr);



  Tcl_CreateCommand(interp, "setCreep", &TclCommand_setCreep, nullptr, nullptr);


  Tcl_CreateCommand(interp, "build", &buildModel, nullptr, nullptr);
  Tcl_CreateCommand(interp, "print", &printModel, nullptr, nullptr);
  Tcl_CreateCommand(interp, "printModel", &printModel, nullptr, nullptr);

  Tcl_CreateCommand(interp, "recorder",          &TclAddRecorder,  domain, nullptr);
  Tcl_CreateCommand(interp, "remove",            &removeObject,    domain, nullptr);

  Tcl_CreateCommand(interp, "findNodeWithID", &findID, domain, nullptr);

// TODO: cmp -- reimplement
//   // Talledo Start
//   Tcl_CreateCommand(interp, "printGID", &printModelGID, nullptr, nullptr);
//   // Talledo End

  Tcl_CreateCommand(interp, "updateElementDomain", &updateElementDomain, nullptr, nullptr);
  Tcl_CreateCommand(interp, "reactions",           &calculateNodalReactions, nullptr, nullptr);
  Tcl_CreateCommand(interp, "nodePressure",        &nodePressure, nullptr, nullptr);
  Tcl_CreateCommand(interp, "nodeBounds",          &nodeBounds, nullptr, nullptr);

  // DAMPING
  Tcl_CreateCommand(interp, "rayleigh",            &rayleighDamping, nullptr, nullptr);
  Tcl_CreateCommand(interp, "setElementRayleighDampingFactors",
                    &setElementRayleighDampingFactors, nullptr, nullptr);

  Tcl_CreateCommand(interp, "region",              &addRegion, nullptr, nullptr);


  Tcl_CreateCommand(interp, "getLoadFactor",       &getLoadFactor, domain, nullptr);
  Tcl_CreateCommand(interp, "localForce",          &localForce,    domain, nullptr);
  Tcl_CreateCommand(interp, "eleType",             &eleType,       domain, nullptr);
  Tcl_CreateCommand(interp, "eleNodes",            &eleNodes,      domain, nullptr);

  Tcl_CreateCommand(interp, "loadConst",           &TclCommand_setLoadConst, domain, nullptr);

  Tcl_CreateCommand(interp, "eleForce",            &eleForce,            domain, nullptr);
  Tcl_CreateCommand(interp, "eleResponse",         &eleResponse,         domain, nullptr);
  Tcl_CreateCommand(interp, "eleDynamicalForce",   &eleDynamicalForce,   domain, nullptr);

  Tcl_CreateCommand(interp, "nodeDOFs",            &nodeDOFs,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeCoord",           &nodeCoord,           domain, nullptr);
  Tcl_CreateCommand(interp, "nodeMass",            &nodeMass,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeVel",             &nodeVel,             domain, nullptr);
  Tcl_CreateCommand(interp, "nodeDisp",            &nodeDisp,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeAccel",           &nodeAccel,           domain, nullptr);
  Tcl_CreateCommand(interp, "nodeResponse",        &nodeResponse,        domain, nullptr);
  Tcl_CreateCommand(interp, "nodeReaction",        &nodeReaction,        domain, nullptr);
  Tcl_CreateCommand(interp, "nodeUnbalance",       &nodeUnbalance,       domain, nullptr);
  Tcl_CreateCommand(interp, "nodeEigenvector",     &nodeEigenvector,     domain, nullptr);

  Tcl_CreateCommand(interp, "setNodeVel",          &setNodeVel,          domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeDisp",         &setNodeDisp,         domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeAccel",        &setNodeAccel,        domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeCoord",        &setNodeCoord,        domain, nullptr);

  Tcl_CreateCommand(interp, "getEleTags",          &getEleTags,          domain, nullptr);
  Tcl_CreateCommand(interp, "getNodeTags",         &getNodeTags,         domain, nullptr);
  Tcl_CreateCommand(interp, "getParamTags",        &getParamTags,        domain, nullptr);
  Tcl_CreateCommand(interp, "getParamValue",       &getParamValue,       domain, nullptr);

  Tcl_CreateCommand(interp, "fixedNodes",          &fixedNodes,          domain, nullptr);
  Tcl_CreateCommand(interp, "fixedDOFs",           &fixedDOFs,           domain, nullptr);
  Tcl_CreateCommand(interp, "constrainedNodes",    &constrainedNodes,    domain, nullptr);
  Tcl_CreateCommand(interp, "constrainedDOFs",     &constrainedDOFs,     domain, nullptr);
  Tcl_CreateCommand(interp, "retainedNodes",       &retainedNodes,       domain, nullptr);
  Tcl_CreateCommand(interp, "retainedDOFs",        &retainedDOFs,        domain, nullptr);

  Tcl_CreateCommand(interp, "getNumElements",      &getNumElements,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleClassTags",     &getEleClassTags,     domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadTags",      &getEleLoadTags,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadData",      &getEleLoadData,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadClassTags", &getEleLoadClassTags, domain, nullptr);


  Tcl_CreateCommand(interp, "sectionForce",        &sectionForce,        domain, nullptr);
  Tcl_CreateCommand(interp, "sectionDeformation",  &sectionDeformation,  domain, nullptr);
  Tcl_CreateCommand(interp, "sectionStiffness",    &sectionStiffness,    domain, nullptr);
  Tcl_CreateCommand(interp, "sectionFlexibility",  &sectionFlexibility,  domain, nullptr);
  Tcl_CreateCommand(interp, "sectionLocation",     &sectionLocation,     domain, nullptr);
  Tcl_CreateCommand(interp, "sectionWeight",       &sectionWeight,       domain, nullptr);
  Tcl_CreateCommand(interp, "basicDeformation",    &basicDeformation,    domain, nullptr);
  Tcl_CreateCommand(interp, "basicForce",          &basicForce,          domain, nullptr);
  Tcl_CreateCommand(interp, "basicStiffness",      &basicStiffness,      domain, nullptr);

  Tcl_CreateCommand(interp, "recorderValue", &OPS_recorderValue, domain, nullptr); // by SAJalali

  // command added for initial state analysis for nDMaterials. Chris McGann, U.Washington
  Tcl_CreateCommand(interp, "InitialStateAnalysis", &InitialStateAnalysis, nullptr, nullptr);


//   TODO: cmp, moved definition to packages/optimization; need to link in optionally
//   Tcl_CreateCommand(interp, "setParameter", &setParameter, nullptr,
//                     nullptr);
  // Tcl_CreateCommand(interp, "sdfResponse",      &sdfResponse, nullptr, nullptr);
  //
  Tcl_CreateCommand(interp, "domainChange", &domainChange, nullptr, nullptr);
  Tcl_CreateCommand(interp, "record",       &TclCommand_record, nullptr, nullptr);
  // Tcl_CreateCommand(interp, "video", &videoPlayer, nullptr, nullptr);
  // Tcl_CreateCommand(interp, "database", &addDatabase, nullptr, nullptr);


  // wipeAnalysis(0, interp, 0, 0);
  return TCL_OK;
}


// by SAJalali
int
OPS_recorderValue(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *domain = G3_getDomain(rt);
  // make sure at least one other argument to contain type of system

  // clmnID starts from 1
  if (argc < 3) {
    opserr << "WARNING want - recorderValue recorderTag clmnID <rowOffset> "
              "<-reset>\n";
    return TCL_ERROR;
  }

  int tag, rowOffset;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> <-reset> "
              "could not read recorderTag \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING recorderValue recorderTag? clmnID - could not read "
              "clmnID \n";
    return TCL_ERROR;
  }
  dof--;
  rowOffset = 0;
  int curArg = 3;
  if (argc > curArg) {
    if (Tcl_GetInt(interp, argv[curArg], &rowOffset) != TCL_OK) {
      opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> "
                "<-reset> could not read rowOffset \n";
      return TCL_ERROR;
    }
    curArg++;
  }
  bool reset = false;
  if (argc > curArg) {
    if (strcmp(argv[curArg], "-reset") == 0)
      reset = true;
    curArg++;
  }
  Recorder *theRecorder = domain->getRecorder(tag);
  double res = theRecorder->getRecordedValue(dof, rowOffset, reset);
  // now we copy the value to the tcl string that is returned
  char buffer[40];
  sprintf(buffer, "%35.8f", res);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
getLoadFactor(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }

  int pattern;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }

  LoadPattern *the_pattern = domain->getLoadPattern(pattern);
  if (the_pattern == 0) {
    opserr << "ERROR load pattern with tag " << pattern
           << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }

  double factor = the_pattern->getLoadFactor();

  char buffer[40];
  sprintf(buffer, "%35.20f", factor);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


// command invoked to build the model, i.e. to invoke buildFE_Model()
// on the ModelBuilder

int
buildModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  ModelBuilder* builder = (ModelBuilder*)G3_getModelBuilder(rt);
  if (!builder)
    builder = theBuilder;
  // TODO: Remove `builtModel` var.
  // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (builder != 0 && builtModel == false) {
    builtModel = true;
    return builder->buildFE_Model();
  } else if (builder != 0 && builtModel == true) {
    opserr << "WARNING Model has already been built - not built again \n";
    return TCL_ERROR;
  } else {
    opserr << "WARNING No ModelBuilder type has been specified \n";
    return TCL_ERROR;
  }
}



int
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char **argv, OPS_Stream &output)
{
  int eleArg = 0;
  if (theAlgorithm == 0)
    return TCL_OK;

  // if just 'print <filename> algorithm'- no flag
  if (argc == 0) {
    theAlgorithm->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> Algorithm flag' get the flag
  int flag;
  if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
    opserr << "WARNING print algorithm failed to get integer flag: \n";
    opserr << argv[eleArg] << endln;
    return TCL_ERROR;
  }
  theAlgorithm->Print(output, flag);
  return TCL_OK;
}



// TODO: consolidate
extern int TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv, EquiSolnAlgo *theAlgorithm);

int
addAlgoRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  if (theAlgorithm != nullptr)
    return TclAddAlgorithmRecorder(clientData, interp, argc, argv, theAlgorithm);
  else
    return TCL_ERROR;
}

/*
int
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      opserr << "WARNING need to specify the commitTag \n";
      return TCL_ERROR;
  }

  if (strcmp(argv[1],"Single") == 0) {
      if (argc < 4) {
        opserr << "WARNING quake single dof motion\n";
        return TCL_ERROR;
      }

      int dof;
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK)
          return TCL_ERROR;

      // read in the ground motion
      GroundMotion *theMotion;
      if (strcmp(argv[3],"ElCentro") == 0) {
          double fact = 1.0;
          if (argc == 5) {
              if (Tcl_GetDouble(interp, argv[4], &fact) != TCL_OK)
                  return TCL_ERROR;
          }
          theMotion = new ElCentroGroundMotion(fact);
      } else {
          opserr << "WARNING quake Single motion - no motion type exists \n";
          return TCL_ERROR;
      }

      Load *theLoad = new SingleExcitation(*theMotion, dof, nextTag++);
      the_domain->addOtherLoad(theLoad);
      return TCL_OK;
  }

  else {
    opserr << "WARNING No quake type exists \n";
    return TCL_ERROR;
  }
}
*/

int
nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeDisp nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeDisp nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeDisp nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, Disp);

  if (nodalResponse == 0)
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned

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
nodeReaction(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeReaction nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeReaction nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeReaction nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, Reaction);

  if (nodalResponse == 0)
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeReaction nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned

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
nodeUnbalance(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeUnbalance nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeUnbalance nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeUnbalance nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, Unbalance);

  if (nodalResponse == 0)
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeUnbalance nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned

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
nodeEigenvector(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int eigenvector = 0;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeEigenvector nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &eigenvector) != TCL_OK) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
    return TCL_ERROR;
  }

  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
      opserr
          << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;
  eigenvector--;
  Node *theNode = domain->getNode(tag);
  const Matrix &theEigenvectors = theNode->getEigenvectors();

  int size = theEigenvectors.noRows();
  int numEigen = theEigenvectors.noCols();

  if (eigenvector < 0 || eigenvector >= numEigen) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
    return TCL_ERROR;
  }

  if (dof >= 0) {
    if (dof >= size) {
      opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = theEigenvectors(dof, eigenvector);

    // now we copy the value to the Tcl string that is returned
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  } else {

    char buffer[40];
    for (int i = 0; i < size; i++) {
      double value = theEigenvectors(i, eigenvector);
      sprintf(buffer, "%35.20f", value);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
eleForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  /*
  Element *theEle = the_domain->getElement(tag);
  if (theEle == 0)
    return TCL_ERROR;

  const Vector &force = theEle->getResistingForce();
  */

  const char *myArgv[1];
  char myArgv0[8];
  strcpy(myArgv0, "forces");
  myArgv[0] = myArgv0;

  const Vector *force = domain->getElementResponse(tag, &myArgv[0], 1);
  if (force != 0) {
    int size = force->Size();

    if (dof >= 0) {

      if (size < dof)
        return TCL_ERROR;

      double value = (*force)(dof);

      // now we copy the value to the tcl string that is returned

      char buffer[40];
      sprintf(buffer, "%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    } else {
      char buffer[40];
      for (int i = 0; i < size; i++) {
        sprintf(buffer, "%35.20f", (*force)(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }
  } else {
    opserr << "WARNING - failed to retrieve element force.\n";
    return TCL_ERROR;
  }
  return TCL_OK;
}

int
localForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - localForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING localForce eleTag? dof? - could not read eleTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING localForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;


  const char *myArgv[1];
  char myArgv0[80];
  strcpy(myArgv0, "localForces");
  myArgv[0] = myArgv0;

  const Vector *force = theDomain->getElementResponse(tag, &myArgv[0], 1);
  if (force != 0) {
    int size = force->Size();

    if (dof >= 0) {

      if (size < dof)
        return TCL_ERROR;

      double value = (*force)(dof);

      // copy the value to the Tcl string that is returned
      char buffer[40];
      sprintf(buffer, "%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    } else {
      char buffer[40];
      for (int i = 0; i < size; i++) {
        sprintf(buffer, "%35.20f", (*force)(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }
  }

  return TCL_OK;
}

int
eleDynamicalForce(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;
  Element *theEle = theDomain->getElement(tag);
  if (theEle == 0)
    return TCL_ERROR;

  const Vector &force = theEle->getResistingForceIncInertia();
  int size = force.Size();

  if (dof >= 0) {

    if (size < dof)
      return TCL_ERROR;

    double value = force(dof);

    // now we copy the value to the tcl string that is returned
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  } else {
    char buffer[40];
    for (int i = 0; i < size; i++) {
      sprintf(buffer, "%35.20f", force(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
eleResponse(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* the_domain = G3_getDomain(rt);
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - eleResponse eleTag? eleArgs...\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  /*
  Element *theEle = the_domain->getElement(tag);
  if (theEle == 0)
    return TCL_ERROR;

  DummyStream dummy;
  Response *theResponse = theEle->setResponse(argv+2, argc-2, dummy);
  if (theResponse == 0) {
    return TCL_ERROR;
  }

  if (theResponse->getResponse() < 0) {
    delete theResponse;
    return TCL_ERROR;
  }

  Information &eleInfo = theResponse->getInformation();
  const Vector &data = eleInfo.getData();
  */

  const Vector *data = the_domain->getElementResponse(tag, argv + 2, argc - 2);
  if (data != 0) {
    int size = data->Size();
    char buffer[40];
    for (int i = 0; i < size; i++) {
      sprintf(buffer, "%f ", (*data)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }
  return TCL_OK;
}

int
findID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - findNodesWithID ?id\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  char buffer[20] = {0};

  while ((theNode = theNodes()) != 0) {
    DOF_Group *theGroup = theNode->getDOF_GroupPtr();
    if (theGroup != 0) {
      const ID &nodeID = theGroup->getID();
      for (int i = 0; i < nodeID.Size(); i++) {
        if (nodeID(i) == tag) {
          sprintf(buffer, "%d ", theNode->getTag());
          Tcl_AppendResult(interp, buffer, NULL);
          break;
        }
      }
    }
  }

  return TCL_OK;
}

int
nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeCoord nodeTag? <dim?>\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeCoord nodeTag? dim? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  int dim = -1;

  if (argc > 2) {
    if (strcmp(argv[2], "X") == 0 || strcmp(argv[2], "x") == 0 ||
        strcmp(argv[2], "1") == 0)
      dim = 0;
    else if (strcmp(argv[2], "Y") == 0 || strcmp(argv[2], "y") == 0 ||
             strcmp(argv[2], "2") == 0)
      dim = 1;
    else if (strcmp(argv[2], "Z") == 0 || strcmp(argv[2], "z") == 0 ||
             strcmp(argv[2], "3") == 0)
      dim = 2;
    else {
      opserr << G3_WARN_PROMPT << "nodeCoord nodeTag? dim? - could not read dim? \n";
      return TCL_ERROR;
    }
  }

  Node *theNode = the_domain->getNode(tag);

  if (theNode == nullptr) {
    opserr << G3_WARN_PROMPT << "Unable to retrieve node with tag '" << tag << "'\n";
    return TCL_ERROR;
  }

  const Vector &coords = theNode->getCrds();

  opserr << "..." << coords;

  char buffer[40];
  int size = coords.Size();
  if (dim == -1) {
    for (int i = 0; i < size; i++) {
      sprintf(buffer, "%35.20f", coords(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
    return TCL_OK;

  } else if (dim < size) {
    double value = coords(dim); // -1 for OpenSees vs C indexing
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  return TCL_ERROR;
}


int
retainedNodes(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;
  bool all = 1;
  int cNode;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
      opserr << "WARNING retainedNodes <cNode?> - could not read cNode? \n";
      return TCL_ERROR;
    }
    all = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = domain->getMPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theMP = mpIter()) != 0) {
    tag = theMP->getNodeRetained();
    if (all || cNode == theMP->getNodeConstrained()) {
      tags.insert(tag);
    }
  }
  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  sort(tagv.begin(), tagv.end());
  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
retainedDOFs(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - retainedDOFs rNode? <cNode?> <cDOF?>\n";
    return TCL_ERROR;
  }

  int rNode;
  if (Tcl_GetInt(interp, argv[1], &rNode) != TCL_OK) {
    opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read "
              "rNode? \n";
    return TCL_ERROR;
  }

  int cNode;
  bool allNodes = 1;
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &cNode) != TCL_OK) {
      opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read "
                "cNode? \n";
      return TCL_ERROR;
    }
    allNodes = 0;
  }

  int cDOF;
  bool allDOFs = 1;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &cDOF) != TCL_OK) {
      opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read "
                "cDOF? \n";
      return TCL_ERROR;
    }
    cDOF--;
    allDOFs = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = domain->getMPs();

  int tag;
  int i;
  int n;
  Vector retained(6);
  while ((theMP = mpIter()) != 0) {
    tag = theMP->getNodeRetained();
    if (tag == rNode) {
      if (allNodes || cNode == theMP->getNodeConstrained()) {
        const ID &rDOFs = theMP->getRetainedDOFs();
        n = rDOFs.Size();
        if (allDOFs) {
          for (i = 0; i < n; i++) {
            retained(rDOFs(i)) = 1;
          }
        } else {
          const ID &cDOFs = theMP->getConstrainedDOFs();
          for (i = 0; i < n; i++) {
            if (cDOF == cDOFs(i))
              retained(rDOFs(i)) = 1;
          }
        }
      }
    }
  }
  char buffer[20];
  for (int i = 0; i < 6; i++) {
    if (retained(i) == 1) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
setNodeCoord(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - setNodeCoord nodeTag? dim? value?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  int dim;
  double value;

  if (Tcl_GetInt(interp, argv[2], &dim) != TCL_OK) {
    opserr
        << "WARNING setNodeCoord nodeTag? dim? value? - could not read dim? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read "
              "value? \n";
    return TCL_ERROR;
  }

  Node *theNode = domain->getNode(tag);

  if (theNode == 0) {
    return TCL_ERROR;
  }

  Vector coords(theNode->getCrds());
  coords(dim - 1) = value;
  theNode->setCrds(coords);

  return TCL_OK;
}

int
updateElementDomain(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char **argv)
{
  // Need to "setDomain" to make the change take effect.
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  ElementIter &theElements = domain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) {
    theElement->setDomain(domain);
  }
  return 0;
}


int
eleType(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - eleType eleTag?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleType eleTag? \n";
    return TCL_ERROR;
  }

  char buffer[80];
  Element *theElement = the_domain->getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING eleType ele " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  const char *type = theElement->getClassType();
  sprintf(buffer, "%s", type);
  Tcl_AppendResult(interp, buffer, NULL);

  return TCL_OK;
}

int
eleNodes(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - eleNodes eleTag?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleNodes eleTag? \n";
    return TCL_ERROR;
  }

  char buffer[20];

  const char *myArgv[1];
  char myArgv0[80];
  strcpy(myArgv0, "nodeTags");
  myArgv[0] = myArgv0;

  // const Vector *tags = the_domain->getElementResponse(tag, &myArgv[0], 1);
  Element *theElement = the_domain->getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING eleNodes ele " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numTags = theElement->getNumExternalNodes();
  const ID &tags = theElement->getExternalNodes();
  for (int i = 0; i < numTags; i++) {
    sprintf(buffer, "%d ", tags(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
nodeDOFs(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[40];
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeDOFs nodeTag?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;
  }


  Node *theNode = the_domain->getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING nodeDOFs node " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();

  DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
  if (theDOFgroup == 0) {
    opserr << "WARNING nodeDOFs DOF group null" << endln;
    return -1;
  }
  const ID &eqnNumbers = theDOFgroup->getID();
  for (int i = 0; i < numDOF; i++) {
    sprintf(buffer, "%d ", eqnNumbers(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
nodeMass(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << "WARNING want - nodeMass nodeTag? nodeDOF?\n";
    return TCL_ERROR;
  }

  int tag, dof;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;
  }

  char buffer[40];

  Node *theNode = the_domain->getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING nodeMass node " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();
  if (dof < 1 || dof > numDOF) {
    opserr << "WARNING nodeMass dof " << dof << " not in range" << endln;
    return TCL_ERROR;
  } else {
    const Matrix &mass = theNode->getMass();
    sprintf(buffer, "%35.20f", mass(dof - 1, dof - 1));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
nodePressure(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING: want - nodePressure nodeTag?\n";
    return TCL_ERROR;
  }
  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING: nodePressure " << argv[1] << "\n";
    return TCL_ERROR;
  }
  double pressure = 0.0;
  Pressure_Constraint *thePC = theDomain.getPressure_Constraint(tag);
  if (thePC != 0) {
    pressure = thePC->getPressure();
  }
  char buffer[80];
  sprintf(buffer, "%35.20f", pressure);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
nodeBounds(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  const int requiredDataSize = 20 * 6;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != 0) {
      delete[] resDataPtr;
    }
    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }

  for (int i = 0; i < requiredDataSize; i++)
    resDataPtr[i] = '\n';

  const Vector &bounds = the_domain->getPhysicalBounds();

  int cnt = 0;
  for (int j = 0; j < 6; j++) {
    cnt += sprintf(&resDataPtr[cnt], "%.6e  ", bounds(j));
  }

  Tcl_SetResult(interp, resDataPtr, TCL_STATIC);

  return TCL_OK;
}

int
nodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeVel nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeVel nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeVel nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = theDomain.getNodeResponse(tag, Vel);

  if (nodalResponse == 0)
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {
    if (size < dof)
      return TCL_ERROR;

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned
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
setNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - setNodeVel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain.getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING setNodeVel -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr
        << "WARNING setNodeVel nodeTag? dof? value?- could not read value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getVel();
    vel(dof) = value;
    theNode->setTrialVel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

int
setNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - setNodeDisp nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain.getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING setNodeDisp -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getDisp();
    vel(dof) = value;
    theNode->setTrialDisp(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}




int
sectionForce(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionForce eleTag? <secNum?> dof? \n";
    return TCL_ERROR;
  }

  int tag, dof;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read "
                "secNum? \n";
      return TCL_ERROR;
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr
        << "WARNING sectionForce eleTag? secNum? dof? - could not read dof? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionForce element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "force";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;
  if (argc < 4) { // For zeroLengthSection
    argcc = 2;
    argvv[1] = c;
  }

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer, "%12.8g", theVec(dof - 1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int
sectionDeformation(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - sectionDeformation eleTag? secNum? dof? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag, secNum, dof;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not "
              "read eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not "
              "read secNum? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not "
              "read dof? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "deformation";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer, "%12.8g", theVec(dof - 1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int
sectionLocation(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionLocation eleTag? secNum? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionLocation eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionLocation eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionLocation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "integrationPoints";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer, "%12.8g", theVec(secNum - 1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int
sectionWeight(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionWeight eleTag? secNum? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING sectionWeight eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr
        << "WARNING sectionWeight eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionWeight element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "integrationWeights";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer, "%12.8g", theVec(secNum - 1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int
sectionStiffness(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionStiffness eleTag? secNum? \n";
    return TCL_ERROR;
  }

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionStiffness eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionStiffness eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionStiffness element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "stiffness";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; i++) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer, "%12.8g ", theMat(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

int
sectionFlexibility(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionFlexibility eleTag? secNum? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionFlexibility eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionFlexibility eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionFlexibility element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "flexibility";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; i++) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer, "%12.8g ", theMat(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

int
basicDeformation(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicDeformation eleTag? \n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum?
  \n"; return TCL_ERROR;
  }
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicDeformation";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int
basicForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicForce eleTag? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicForce eleTag? dofNum? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum?
  \n"; return TCL_ERROR;
  }
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicForce";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int
basicStiffness(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicStiffness eleTag? \n";
    return TCL_ERROR;
  }

  // opserr << "sectionDeformation: ";
  // for (int i = 0; i < argc; i++)
  //  opserr << argv[i] << ' ' ;
  // opserr << endln;

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicStiffness eleTag? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum?
  \n"; return TCL_ERROR;
  }
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicStiffness element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicStiffness";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMatrix = *(info.theMatrix);
  int nbf = theMatrix.noCols();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nbf; j++) {
      sprintf(buffer, "%12.8f ", theMatrix(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

// added by C.McGann, U.Washington
int
InitialStateAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char **argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING: Incorrect number of arguments for InitialStateAnalysis "
              "command"
           << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "on") == 0) {
    opserr << "InitialStateAnalysis ON" << endln;

    // set global variable to true
    // FMK changes for parallel:
    // ops_InitialStateAnalysis = true;

    Parameter *theP = new InitialStateParameter(true);
    theDomain.addParameter(theP);
    delete theP;

    return TCL_OK;

  } else if (strcmp(argv[1], "off") == 0) {
    opserr << "InitialStateAnalysis OFF" << endln;

    // call revert to start to zero the displacements
    theDomain.revertToStart();

    // set global variable to false
    // FMK changes for parallel
    // ops_InitialStateAnalysis = false;
    Parameter *theP = new InitialStateParameter(false);
    theDomain.addParameter(theP);
    delete theP;

    return TCL_OK;

  } else {
    opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or "
              "InitialStateAnalysis off"
           << endln;

    return TCL_ERROR;
  }
}

int
rayleighDamping(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  if (argc < 5) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - not enough "
              "arguments to command\n";
    return TCL_ERROR;
  }
  double alphaM, betaK, betaK0, betaKc;
  if (Tcl_GetDouble(interp, argv[1], &alphaM) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read alphaM? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &betaK) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &betaK0) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK0? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &betaKc) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaKc? \n";
    return TCL_ERROR;
  }

  Domain *the_domain = G3_getDomain(G3_getRuntime(interp));
  the_domain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

  return TCL_OK;
}

int
setElementRayleighDampingFactors(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv)
{
  if (argc < 6) {
    opserr << "WARNING setElementRayleighDampingFactors eleTag? alphaM? betaK? "
              "betaK0? betaKc? - not enough arguments to command\n";
    return TCL_ERROR;
  }
  int eleTag;
  double alphaM, betaK, betaK0, betaKc;

  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read eleTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read alphaM? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &betaK0) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK0? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &betaKc) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaKc? \n";
    return TCL_ERROR;
  }

  Element *theEle = theDomain.getElement(eleTag);
  theEle->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}

extern int TclAddMeshRegion(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char **argv, Domain &theDomain);

int
addRegion(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Domain *the_domain = G3_getDomain(G3_getRuntime(interp));
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, the_domain);
  return TclAddMeshRegion(clientData, interp, argc, argv, theDomain);
}




int
getNumElements(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char **argv)
{
  char buffer[20];

  sprintf(buffer, "%d ", theDomain.getNumElements());
  Tcl_AppendResult(interp, buffer, NULL);

  return TCL_OK;
}

int
getEleClassTags(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{

  if (argc == 1) {
    Element *theEle;
    ElementIter &eleIter = theDomain.getElements();

    char buffer[20];

    while ((theEle = eleIter()) != 0) {
      sprintf(buffer, "%d ", theEle->getClassTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }
  } else if (argc == 2) {
    int eleTag;

    if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
      opserr << "WARNING getParamValue -- could not read paramTag \n";
      return TCL_ERROR;
    }

    Element *theEle = theDomain.getElement(eleTag);

    char buffer[20];

    sprintf(buffer, "%d ", theEle->getClassTag());
    Tcl_AppendResult(interp, buffer, NULL);

  } else {
    opserr << "WARNING want - getEleClassTags <eleTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadClassTags(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char **argv)
{

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = theDomain.getLoadPatterns();

    char buffer[20];

    while ((thePattern = thePatterns()) != 0) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != 0) {
        sprintf(buffer, "%d ", theLoad->getClassTag());
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << "WARNING getEleLoadClassTags -- could not read patternTag\n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = theDomain.getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << "ERROR load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadClassTags\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    char buffer[20];

    while ((theLoad = theEleLoads()) != 0) {
      sprintf(buffer, "%d ", theLoad->getClassTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << "WARNING want - getEleLoadClassTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadTags(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char **argv)
{

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = theDomain.getLoadPatterns();

    char buffer[20];

    while ((thePattern = thePatterns()) != 0) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != 0) {
        sprintf(buffer, "%d ", theLoad->getElementTag());
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << "WARNING getEleLoadTags -- could not read patternTag \n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = theDomain.getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << "ERROR load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadTags\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    char buffer[20];

    while ((theLoad = theEleLoads()) != 0) {
      sprintf(buffer, "%d ", theLoad->getElementTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << "WARNING want - getEleLoadTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadData(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char **argv)
{

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = theDomain.getLoadPatterns();

    char buffer[40];
    int typeEL;

    while ((thePattern = thePatterns()) != 0) {
      ElementalLoadIter &theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != 0) {
        const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

        int eleLoadDataSize = eleLoadData.Size();
        opserr << "eleLoadDataSize: " << eleLoadDataSize << "\n";
        for (int i = 0; i < eleLoadDataSize; i++) {
          sprintf(buffer, "%35.20f ", eleLoadData(i));
          Tcl_AppendResult(interp, buffer, NULL);
        }
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << "WARNING getEleLoadData -- could not read patternTag \n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = theDomain.getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << "ERROR load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadData\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    int typeEL;
    char buffer[40];

    while ((theLoad = theEleLoads()) != 0) {
      const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

      int eleLoadDataSize = eleLoadData.Size();
      for (int i = 0; i < eleLoadDataSize; i++) {
        sprintf(buffer, "%35.20f ", eleLoadData(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else {
    opserr << "WARNING want - getEleLoadTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Element *theEle;
  ElementIter &eleIter = theDomain.getElements();

  char buffer[20];

  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getNodeTags(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *the_domain = G3_getDomain(rt);
  Node *node;
  if (the_domain==nullptr)
    return TCL_ERROR;

  NodeIter &nodeIter = the_domain->getNodes();

  char buffer[20];

  while ((node = nodeIter()) != 0) {
    sprintf(buffer, "%d ", node->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getParamTags(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{
  Parameter *theEle;
  ParameterIter &eleIter = theDomain.getParameters();

  char buffer[20];

  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getParamValue(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "Insufficient arguments to getParamValue" << endln;
    return TCL_ERROR;
  }

  int paramTag;

  if (Tcl_GetInt(interp, argv[1], &paramTag) != TCL_OK) {
    opserr << "WARNING getParamValue -- could not read paramTag \n";
    return TCL_ERROR;
  }

  Parameter *theEle = theDomain.getParameter(paramTag);

  char buffer[40];

  sprintf(buffer, "%35.20f", theEle->getValue());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
TclCommand_record(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);
  ((Domain*)clientData)->record(false);
  return TCL_OK;
}

int
elementActivate(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  int eleTag;
  int argLoc = 1;
  int Nelements = argc;
  ID activate_us(0, Nelements);

  while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
    activate_us.insert(eleTag);
    ++argLoc;
  }

  theDomain.activateElements(activate_us);

  return TCL_OK;
}

int
elementDeactivate(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv)
{

  int eleTag;
  int argLoc = 1;
  int Nelements = argc;
  ID deactivate_us(0, Nelements);

  while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
    deactivate_us.insert(eleTag);
    ++argLoc;
  }

  theDomain.deactivateElements(deactivate_us);
  return TCL_OK;
}
