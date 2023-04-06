//
// Written: fmk
// Created: 07/99
//
// Description: This file contains the class definition for TclBasicBuilder.
// A TclBasicBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework. currently these elements include:
//
#include <stdlib.h>
#include <string.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>

#include <RigidRod.h>
#include <RigidBeam.h>

#include <CrdTransf.h>

#include <NodalLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam2dPartialUniformLoad.h>
#include <Beam2dTempLoad.h>
#include <Beam2dThermalAction.h>  //L.Jiang [SIF]
#include <Beam3dThermalAction.h>  //L.Jiang [SIF]
#include <ShellThermalAction.h>   //L.Jiang [SIF]
#include <ThermalActionWrapper.h> //L.Jiang [SIF]
#include <NodalThermalAction.h>   //L.Jiang [SIF]

#include <Beam3dPointLoad.h>
#include <Beam3dUniformLoad.h>
#include <Beam3dPartialUniformLoad.h>
#include <BrickSelfWeight.h>
#include <SurfaceLoader.h>
#include <SelfWeight.h>
#include <LoadPattern.h>

#include <SectionForceDeformation.h>
#include <SectionRepres.h>

#include <NDMaterial.h>
#include <TclBasicBuilder.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>

#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>                   //L.Jiang [SIF]
#include <vector>                                    //L.Jiang [SIF]
using std::vector;                                   // L.Jiang [SIF]

// Added by Scott J. Brandenberg (sjbrandenberg@ucdavis.edu)
#include <PySimple1Gen.h>
#include <TzSimple1Gen.h>
// End added by SJB

// Added by Prishati Raychowdhury  (PRC)
#include <ShallowFoundationGen.h>
// end PRC

#include <YieldSurface_BC.h>
#include <YS_Evolution.h>
#include <PlasticHardeningMaterial.h>

#ifdef OPSDEF_DAMAGE
#  include <DamageModel.h> //!!
#endif

#include <FrictionModel.h>

#include <StiffnessDegradation.h>
#include <UnloadingRule.h>
#include <StrengthDegradation.h>
// #include <HystereticBackbone.h>
#include <BeamIntegration.h>

////////////////////// gnp adding damping
#include <Element.h>
////////////////////////////////////////////

extern int OPS_ResetInput(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *domain,
                          TclBuilder *builder);
#include <packages.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Domain *theTclDomain = 0;
static TclBasicBuilder *theTclBuilder = 0;

extern LoadPattern *theTclLoadPattern;
// extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
static int eleLoadTag = 0;

//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//
// REMO
extern int TclCommand_addPatch(ClientData, Tcl_Interp*,
                               int argc, TCL_Char **const);

extern int TclCommand_addFiber(ClientData, Tcl_Interp*,
                               int argc, TCL_Char **const);

extern int TclCommand_addReinfLayer(ClientData, Tcl_Interp*,
                                   int argc, TCL_Char **const);




static int TclCommand_addParameter(ClientData, Tcl_Interp*,
                                   int argc, TCL_Char ** const);

extern int TclCommand_addElement(ClientData, Tcl_Interp*,
                                 int argc, TCL_Char ** const);

static int TclCommand_mesh(ClientData, Tcl_Interp*, int argc,
                           TCL_Char ** const);

static int TclCommand_remesh(ClientData, Tcl_Interp*,
                             int argc, TCL_Char ** const);


static int TclCommand_addBeamIntegration(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);

static int TclCommand_addLimitCurve(ClientData, Tcl_Interp*,
                                    int argc, TCL_Char ** const);

extern int TclCommand_addNDMaterial(ClientData, Tcl_Interp*,
                                    int argc, TCL_Char ** const);


extern int TclCommand_addSection(ClientData, Tcl_Interp*,
                                 int argc, TCL_Char ** const);

static int TclCommand_addYieldSurface_BC(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);

static int TclCommand_addYS_EvolutionModel(ClientData,
                                           Tcl_Interp*, int argc,
                                           TCL_Char ** const);

static int TclCommand_addYS_PlasticMaterial(ClientData,
                                            Tcl_Interp*, int argc,
                                            TCL_Char ** const);

int
TclCommand_addCyclicModel(ClientData, Tcl_Interp*, int argc,
                          TCL_Char ** const);

#ifdef OPSDEF_DAMAGE
int
TclCommand_addDamageModel(ClientData, Tcl_Interp*, int argc,
                          TCL_Char ** const);
#endif // OPSDEF_DAMAGE
/*
static int TclCommand_addTimeSeries(ClientData, Tcl_Interp*,
                                    int argc, TCL_Char ** const);

static int TclCommand_addPattern(ClientData, Tcl_Interp*,
                                 int argc, TCL_Char ** const);

static int TclCommand_addSeries(ClientData, Tcl_Interp*,
                                int argc, TCL_Char ** const);
*/
static int TclCommand_addHomogeneousBC(ClientData,
                                       Tcl_Interp*, int argc,
                                       TCL_Char ** const);
static int TclCommand_addHomogeneousBC_X(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);
static int TclCommand_addHomogeneousBC_Y(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);
static int TclCommand_addHomogeneousBC_Z(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);
static int TclCommand_addEqualDOF_MP(ClientData, Tcl_Interp*,
                                     int argc, TCL_Char ** const);

static int TclCommand_addEqualDOF_MP_Mixed(ClientData,
                                           Tcl_Interp*, int argc,
                                           TCL_Char ** const);

static int TclCommand_RigidLink(ClientData, Tcl_Interp*,
                                int argc, TCL_Char ** const);

static int TclCommand_addMP(ClientData, Tcl_Interp*, int argc,
                            TCL_Char ** const);

static int TclCommand_addNodalLoad(ClientData, Tcl_Interp*,
                                   int argc, TCL_Char ** const);


static int TclCommand_addNodalMass(ClientData, Tcl_Interp*,
                                   int argc, TCL_Char ** const);
static int TclCommand_addSP(ClientData, Tcl_Interp*, int argc,
                            TCL_Char ** const);

static int TclCommand_addImposedMotionSP(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);
// Added by Scott J. Brandenberg
static int TclCommand_doPySimple1Gen(ClientData, Tcl_Interp*,
                                     int argc, TCL_Char ** const);

static int TclCommand_doTzSimple1Gen(ClientData, Tcl_Interp*,
                                     int argc, TCL_Char ** const);
// End added by SJB

// Added by Prishati Raychowdhury (UCSD)
static int TclBasicBuilder_doShallowFoundationGen(ClientData,
                                                  Tcl_Interp*, int argc,
                                                  TCL_Char ** const);
// End PRC

// Leo
static int TclBasicBuilder_addRemoHFiber(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);

static int TclCommand_addFrictionModel(ClientData,
                                       Tcl_Interp*, int argc,
                                       TCL_Char ** const);

static int TclCommand_addStiffnessDegradation(ClientData,
                                              Tcl_Interp*, int argc,
                                              TCL_Char ** const);

static int TclCommand_addUnloadingRule(ClientData,
                                       Tcl_Interp*, int argc,
                                       TCL_Char ** const);

static int TclCommand_addStrengthDegradation(ClientData,
                                             Tcl_Interp*, int argc,
                                             TCL_Char ** const);
/*
static int TclCommand_addHystereticBackbone(ClientData,
                                            Tcl_Interp*, int argc,
                                            TCL_Char ** const);
*/

extern int TclCommand_addGroundMotion(ClientData, Tcl_Interp*,
                                      int argc, TCL_Char ** const);

/// added by ZHY
static int TclCommand_UpdateMaterialStage(ClientData,
                                          Tcl_Interp*, int argc,
                                          TCL_Char ** const);

static Tcl_CmdProc TclCommand_UpdateMaterials;

extern Tcl_CmdProc TclBasicBuilderUpdateParameterCommand;

/// added by ZHY

////////////////gnp adding rayleigh //////////////////////////
static int TclCommand_addElementRayleigh(ClientData,
                                         Tcl_Interp*, int argc,
                                         TCL_Char ** const);
///////////////////////////////////////////////////////////////


extern int TclCommand_addHFiber(ClientData, Tcl_Interp*,
                                int argc, TCL_Char **const,
                                TclBasicBuilder *theTclBuilder);
/*
extern int TclCommand_addReinfLayer(ClientData, Tcl_Interp*,
                                    int argc, TCL_Char **,
                                    TclBasicBuilder *theTclBuilder);
*/
extern int TclCommand_addGeomTransf(ClientData, Tcl_Interp *, int, TCL_Char ** const);

static int TclCommand_Package(ClientData, Tcl_Interp*,
                              int argc, TCL_Char ** const);


// Added by Alborz Ghofrani - U.Washington
extern int TclCommand_GenerateInterfacePoints(ClientData,
                                              Tcl_Interp*, int argc,
                                              TCL_Char ** const);
// End Added by Alborz


//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclBasicBuilder::TclBasicBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM,
                                 int NDF)
    : TclBuilder(theDomain, NDM, NDF), theInterp(interp)
{
  theSections = new ArrayOfTaggedObjects(32);
  theSectionRepresents = new ArrayOfTaggedObjects(32);

  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "parameter", TclCommand_addParameter, NULL, NULL);
  Tcl_CreateCommand(interp, "addToParameter", TclCommand_addParameter, NULL, NULL);
  Tcl_CreateCommand(interp, "updateParameter", TclCommand_addParameter, NULL, NULL);

  Tcl_CreateCommand(interp, "mesh", TclCommand_mesh, NULL, NULL);
  Tcl_CreateCommand(interp, "remesh", TclCommand_remesh, NULL, NULL);


  Tcl_CreateCommand(interp, "load", TclCommand_addNodalLoad, NULL, NULL);

  Tcl_CreateCommand(interp, "imposedMotion", TclCommand_addImposedMotionSP, NULL, NULL);

  Tcl_CreateCommand(interp, "imposedSupportMotion",
                    TclCommand_addImposedMotionSP, NULL, NULL);

  Tcl_CreateCommand(interp, "groundMotion", TclCommand_addGroundMotion, NULL, NULL);

  Tcl_CreateCommand(interp, "equalDOF",       TclCommand_addEqualDOF_MP, NULL, NULL);
  Tcl_CreateCommand(interp, "equalDOF_Mixed", TclCommand_addEqualDOF_MP_Mixed, NULL, NULL);
  Tcl_CreateCommand(interp, "rigidLink",      &TclCommand_RigidLink, NULL, NULL);

  Tcl_CreateCommand(interp, "sp", TclCommand_addSP, NULL, NULL);
  Tcl_CreateCommand(interp, "mp", TclCommand_addMP, NULL, NULL);

  Tcl_CreateCommand(interp, "PySimple1Gen", TclCommand_doPySimple1Gen, NULL, NULL);
  Tcl_CreateCommand(interp, "TzSimple1Gen", TclCommand_doTzSimple1Gen, NULL, NULL);

  Tcl_CreateCommand(interp, "ShallowFoundationGen",
                    TclBasicBuilder_doShallowFoundationGen, NULL, NULL);

  // Added by LEO
  Tcl_CreateCommand(interp, "Hfiber",               TclBasicBuilder_addRemoHFiber, NULL, NULL);

  Tcl_CreateCommand(interp, "frictionModel",        TclCommand_addFrictionModel, NULL, NULL);
  Tcl_CreateCommand(interp, "beamIntegration",  TclCommand_addBeamIntegration,  NULL, NULL);
  Tcl_CreateCommand(interp, "yieldSurface_BC", TclCommand_addYieldSurface_BC, NULL, NULL);
  Tcl_CreateCommand(interp, "ysEvolutionModel", TclCommand_addYS_EvolutionModel, NULL, NULL);
  Tcl_CreateCommand(interp, "plasticMaterial", TclCommand_addYS_PlasticMaterial, NULL, NULL);
  Tcl_CreateCommand(interp, "cyclicModel", TclCommand_addCyclicModel, NULL, NULL);
  Tcl_CreateCommand(interp, "stiffnessDegradation", TclCommand_addStiffnessDegradation, NULL, NULL);
  Tcl_CreateCommand(interp, "unloadingRule",        TclCommand_addUnloadingRule, NULL, NULL);
  Tcl_CreateCommand(interp, "strengthDegradation",  TclCommand_addStrengthDegradation, NULL, NULL);
  // Tcl_CreateCommand(interp, "updateMaterialStage",  TclCommand_UpdateMaterialStage, NULL, NULL);
  // Tcl_CreateCommand(interp, "updateMaterials",      TclCommand_UpdateMaterials, NULL, NULL);
  Tcl_CreateCommand(interp, "loadPackage",          TclCommand_Package, NULL, NULL);
  Tcl_CreateCommand(interp, "setElementRayleighFactors", TclCommand_addElementRayleigh, NULL, NULL);


  // set the static pointers in this file
  theTclBuilder = this;
  theTclDomain = &theDomain;
  theTclLoadPattern = 0;
  // theTclMultiSupportPattern = 0;

  nodeLoadTag = 0;
  eleArgStart = 0;
  Tcl_SetAssocData(interp, "OPS::theTclBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theTclDomain", NULL, (ClientData)&theDomain);
}

TclBasicBuilder::~TclBasicBuilder()
{
  theSections->clearAll();
  theSectionRepresents->clearAll();
  delete theSections;
  delete theSectionRepresents;

  // set the pointers to 0
  theTclDomain = nullptr;
  theTclBuilder = nullptr;
  theTclLoadPattern = nullptr;
  // theTclMultiSupportPattern = 0;

  // may possibly invoke Tcl_DeleteCommand() later
  Tcl_DeleteCommand(theInterp, "parameter");
  Tcl_DeleteCommand(theInterp, "addToParameter");
  Tcl_DeleteCommand(theInterp, "updateParameter");
  Tcl_DeleteCommand(theInterp, "mesh");
  Tcl_DeleteCommand(theInterp, "remesh");
  Tcl_DeleteCommand(theInterp, "background");
  Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  Tcl_DeleteCommand(theInterp, "imposedSupportMotion");
  Tcl_DeleteCommand(theInterp, "groundMotion");
  Tcl_DeleteCommand(theInterp, "equalDOF");
  Tcl_DeleteCommand(theInterp, "sp");
  Tcl_DeleteCommand(theInterp, "mp");
  Tcl_DeleteCommand(theInterp, "PySimple1Gen"); // Added by Scott J. Brandenberg
  Tcl_DeleteCommand(theInterp, "TzSimple1Gen"); // Added by Scott J. Brandenberg

  Tcl_DeleteCommand(theInterp, "fiber");
  Tcl_DeleteCommand(theInterp, "Hfiber"); // LEO
  Tcl_DeleteCommand(theInterp, "updateMaterialStage");
  Tcl_DeleteCommand(theInterp, "updateMaterials");

  Tcl_DeleteCommand(theInterp, "frictionModel");
  Tcl_DeleteCommand(theInterp, "unloadingRule");
  Tcl_DeleteCommand(theInterp, "stiffnessDegradation");
  Tcl_DeleteCommand(theInterp, "strengthDegradation");
  Tcl_DeleteCommand(theInterp, "hystereticBackbone");

  Tcl_DeleteCommand(theInterp, "damageModel");

  Tcl_DeleteCommand(theInterp, "loadPackage");
  Tcl_DeleteCommand(
      theInterp,
      "generateInterfacePoints"); // Added by Alborz Ghofrani - U.Washington
}

//
// CLASS METHODS
//

int
TclBasicBuilder::addSection(SectionForceDeformation &theSection)
{
  //  bool result = theSections->addComponent(&theSection);
  bool result = OPS_addSectionForceDeformation(&theSection);
  if (result == true)
    return 0;
  else {
    opserr << "TclBasicBuilder::addSection() - failed to add section: "
           << theSection;
    return -1;
  }
}

SectionForceDeformation *
TclBasicBuilder::getSection(int tag)
{
  return OPS_getSectionForceDeformation(tag);
}


int
TclBasicBuilder::addSectionRepres(SectionRepres &theSectionRepres)
{
  bool result = theSectionRepresents->addComponent(&theSectionRepres);

  if (result == true)
    return 0;
  else {
    opserr << "TclBasicBuilder::addSectionRepres() - failed to add "
              "SectionRepres\n";
    return -1;
  }
}

SectionRepres *
TclBasicBuilder::getSectionRepres(int tag)
{
  TaggedObject *mc = theSectionRepresents->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  SectionRepres *result = (SectionRepres *)mc;
  return result;
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclCommand_addElementRayleigh(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{

  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  // make sure corect number of arguments on command line
  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: setElementRayleighFactors elementTag?  alphaM? $betaK? "
              "$betaKinit? $betaKcomm? \n";
    return TCL_ERROR;
  }

  int eleTag = 0;

  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << "WARNING: setElementRayleighFactors invalid eleTag: " << argv[1];
    opserr << " \n";
    return TCL_ERROR;
  }

  double alphaM, betaK, betaKinit, betaKcomm;

  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "alphaM: " << argv[2] << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaK: " << argv[3] << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &betaKinit) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKinit: " << argv[4] << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &betaKcomm) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKcomm: " << argv[5] << endln;
    return TCL_ERROR;
  }

  Element *elePtr = theTclDomain->getElement(eleTag);

  if (elePtr == 0)
    opserr << "WARNING : setElementRayleighFactors invalid eleTag: " << eleTag
           << " the element does not exist in the domain \n";

  if (elePtr->setRayleighDampingFactors(alphaM, betaK, betaKinit, betaKcomm) !=
      0) {
    opserr << "ERROR : setElementRayleighFactors: FAILED to add damping "
              "factors for element "
           << eleTag << "\n";
  }

  return TCL_OK;
}
/////////////////////////////   gnp adding element damping

// the function for creating ne material objects and patterns is in a seperate
// file. this allows new material and patternobjects to be added without
// touching this file. does so at the expense of an extra procedure call.

extern int TclBasicBuilderParameterCommand(ClientData clientData,
                                           Tcl_Interp *interp, int argc,
                                           TCL_Char ** const argv, Domain *theDomain,
                                           TclBasicBuilder *theTclBuilder);
int
TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)

{
  return TclBasicBuilderParameterCommand(clientData, interp, argc, argv,
                                         theTclDomain, theTclBuilder);
}


// extern int OPS_LineMesh(Domain& domain, int ndm);
// extern int OPS_TriMesh(Domain& domain);
// extern int OPS_TriReMesh(Domain& domain, int ndf);
int
TclCommand_mesh(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  // make sure corect number of arguments on command line
  if (argc < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: mesh type? ...>\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
                 theTclBuilder);

  // mesh type
  int res = 0;
  // if (strcmp(argv[1], "line") == 0) {
  // 	res = OPS_LineMesh(*theTclDomain,ndm);
  // } else if (strcmp(argv[1], "tri") == 0) {
  // 	res = OPS_TriMesh(*theTclDomain);
  // } else {
  // 	opserr<<"WARNING: mesh type "<<argv[1]<<" is unknown\n";
  // 	return TCL_ERROR;
  // }

  if (res < 0) {
    return TCL_ERROR;
  }

  return 0;
}

int
TclCommand_remesh(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == nullptr) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }


  // make sure corect number of arguments on command line
  if (argc < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: mesh type? ...>\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
                 theTclBuilder);

  // mesh type
  int res = 0;
  // if (strcmp(argv[1], "line") == 0) {
  // 	//res = OPS_LineMesh(*theTclDomain,ndm);
  // } else if (strcmp(argv[1], "tri") == 0) {
  // 	res = OPS_TriReMesh(*theTclDomain,ndf);
  // } else {
  // 	opserr<<"WARNING: remesh type "<<argv[1]<<" is unknown\n";
  // 	return TCL_ERROR;
  // }

  if (res < 0) {
    return TCL_ERROR;
  }

  return 0;
}

#if defined(OPSDEF_Element_PFEM)
extern int OPS_BgMesh();

int
TclCommand_backgroundMesh(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 1, argc, argv, theTclDomain,
                 theTclBuilder);

  if (OPS_BgMesh() >= 0)
    return TCL_OK;
  else
    return TCL_ERROR;
  return TCL_OK;
}
#endif // _OPS_Element_PFEM

extern void *OPS_LobattoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LegendreBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_NewtonCotesBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_RadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_TrapezoidalBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_CompositeSimpsonBeamIntegration(int &integrationTag,
                                                 ID &secTags);
extern void *OPS_UserDefinedBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_FixedLocationBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LowOrderBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_MidDistanceBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_UserHingeBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeMidpointBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauTwoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeEndpointBeamIntegration(int &integrationTag, ID &secTags);

int
TclCommand_addBeamIntegration(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  if (argc < 2) {
    opserr << "WARNING: want beamIntegration type itag...\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
                 theTclBuilder);

  int iTag;
  ID secTags;
  BeamIntegration *bi = 0;
  if (strcmp(argv[1], "Lobatto") == 0) {
    bi = (BeamIntegration *)OPS_LobattoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Legendre") == 0) {
    bi = (BeamIntegration *)OPS_LegendreBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "NewtoCotes") == 0) {
    bi = (BeamIntegration *)OPS_NewtonCotesBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Radau") == 0) {
    bi = (BeamIntegration *)OPS_RadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Trapezoidal") == 0) {
    bi = (BeamIntegration *)OPS_TrapezoidalBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "CompositeSimpson") == 0) {
    bi = (BeamIntegration *)OPS_CompositeSimpsonBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserDefined") == 0) {
    bi = (BeamIntegration *)OPS_UserDefinedBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "FixedLocation") == 0) {
    bi = (BeamIntegration *)OPS_FixedLocationBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "LowOrder") == 0) {
    bi = (BeamIntegration *)OPS_LowOrderBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "MidDistance") == 0) {
    bi = (BeamIntegration *)OPS_MidDistanceBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserHinge") == 0) {
    bi = (BeamIntegration *)OPS_UserHingeBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeMidpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeMidpointBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadau") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadauTwo") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauTwoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeEndpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeEndpointBeamIntegration(iTag, secTags);
  } else {
    opserr << "WARNING: integration type " << argv[1] << " is unknown\n";
    return TCL_ERROR;
  }

  if (bi == 0) {
    opserr << "WARNING: failed to create beam integration\n";
    return TCL_ERROR;
  }

  BeamIntegrationRule *rule = new BeamIntegrationRule(iTag, bi, secTags);
  if (rule == 0) {
    opserr << "WARNING: failed to create beam integration\n";
    delete bi;
    return TCL_ERROR;
  }

  // Now add the
  if (OPS_addBeamIntegrationRule(rule) == false) {
    opserr << "WARNING: could not add BeamIntegrationRule.";
    delete rule; // invoke the destructor, otherwise mem leak
    return TCL_ERROR;
    ;
  }

  return TCL_OK;
}

// extern int Tcl_AddLimitCurveCommand(ClientData clienData, Tcl_Interp *interp,
//                                     int argc, TCL_Char ** const argv,
//                                     Domain *theDomain);
// 
// int
// TclCommand_addLimitCurve(ClientData clientData, Tcl_Interp *interp, int argc,
//                          TCL_Char ** const argv)
// {
//   return Tcl_AddLimitCurveCommand(clientData, interp, argc, argv, theTclDomain);
// }
//
extern int
TclBasicBuilderYieldSurface_BCCommand(ClientData clienData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv,
                                      TclBasicBuilder *theTclBuilder);

int
TclCommand_addYieldSurface_BC(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)

{
  return TclBasicBuilderYieldSurface_BCCommand(clientData, interp, argc, argv,
                                               theTclBuilder);
}

extern int TclBasicBuilderYS_EvolutionModelCommand(
    ClientData clienData, Tcl_Interp *interp, int argc, TCL_Char ** const argv,
    TclBasicBuilder *theTclBuilder);

int
TclCommand_addYS_EvolutionModel(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)

{
  return TclBasicBuilderYS_EvolutionModelCommand(clientData, interp, argc, argv,
                                                 theTclBuilder);
}

extern int
TclBasicBuilderPlasticMaterialCommand(ClientData clienData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv,
                                      TclBasicBuilder *theTclBuilder);

int
TclCommand_addYS_PlasticMaterial(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)

{
  return TclBasicBuilderPlasticMaterialCommand(clientData, interp, argc, argv,
                                               theTclBuilder);
}

//!!
extern int TclBasicBuilderCyclicModelCommand(ClientData clienData,
                                             Tcl_Interp *interp, int argc,
                                             TCL_Char ** const argv,
                                             TclBasicBuilder *theTclBuilder);
int
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)

{
  return TclBasicBuilderCyclicModelCommand(clientData, interp, argc, argv,
                                           theTclBuilder);
}

extern int TclBasicBuilderDamageModelCommand(ClientData clienData,
                                             Tcl_Interp *interp, int argc,
                                             TCL_Char ** const argv);

#ifdef OPSDEF_DAMAGE
int
TclCommand_addDamageModel(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)

{
  return TclBasicBuilderDamageModelCommand(clientData, interp, argc, argv);
}
#endif // OPSDEF_DAMAGE


int
TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";
    return TCL_ERROR;
  }

  //  int ndf = builder->getNDF();

  int ndf = argc - 2;
  NodalLoad *theLoad = 0;

  bool isLoadConst = false;
  bool userSpecifiedPattern = false;
  int loadPatternTag = 0;
  // The above definition are moved forward for the use in both cases

  //-------------Adding Proc for NodalThermalAction, By Liming Jiang, [SIF] 2017
  if ((strcmp(argv[2], "-NodalThermal") == 0) ||
      (strcmp(argv[2], "-nodalThermal") == 0)) {

#if 0
    int nodeId;
    if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      opserr << "WARNING invalid nodeId: " << argv[1] << endln;
      return TCL_ERROR;
    }

    Vector *thecrds = new Vector();
    Node *theNode = theTclDomain->getNode(nodeId);
    if (theNode == 0) {
      opserr << "WARNING invalid nodeID: " << argv[1] << endln;
      return TCL_ERROR;
    }
    (*thecrds) = theNode->getCrds();

    int count = 3;
    if (strcmp(argv[count], "-source") == 0) {
      count++;
      const char *pwd = getInterpPWD(interp);
      simulationInfo.addInputFile(argv[count], pwd);
      TimeSeries *theSeries;

      int dataLen =
          9; // default num of temperature input for nodal ThermalAction;

      if (argc - count == 5) {
        // which indicates the nodal thermal action is applied to 3D I section
        // Beam;
        dataLen = 15;
        theSeries = new PathTimeSeriesThermal(nodeId, argv[count], dataLen);
        count++;
        double RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4;
        if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
          opserr << "WARNING NodalLoad - invalid loc1  " << argv[count]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
          opserr << "WARNING NodalLoad - invalid loc2  " << argv[count + 1]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 2], &RcvLoc3) != TCL_OK) {
          opserr << "WARNING NodalLoad - invalid loc3  " << argv[count + 2]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 3], &RcvLoc4) != TCL_OK) {
          opserr << "WARNING NodalLoad - invalid loc4  " << argv[count + 3]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        // end of recieving data;
        theLoad = new NodalThermalAction(nodeLoadTag, nodeId, RcvLoc1, RcvLoc2,
                                         RcvLoc3, RcvLoc4, theSeries, thecrds);
      }
      // end of for 15 data input;
      else if (argc - count == 3 || argc - count == 10) {

        theSeries = new PathTimeSeriesThermal(nodeId, argv[count]);
        count++;
        Vector locy;
        if (argc - count == 2) {
          double RcvLoc1, RcvLoc2;
          if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
            opserr << "WARNING NodalLoad - invalid loc1  " << argv[count]
                   << " for NodalThermalAction\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
            opserr << "WARNING NodalLoad - invalid loc2  " << argv[count + 1]
                   << " for NodalThermalAction\n";
            return TCL_ERROR;
          }
          locy = Vector(9);
          locy(0) = RcvLoc1;
          locy(1) = (7 * RcvLoc1 + 1 * RcvLoc2) / 8;
          locy(2) = (6 * RcvLoc1 + 2 * RcvLoc2) / 8;
          locy(3) = (5 * RcvLoc1 + 3 * RcvLoc2) / 8;
          locy(4) = (4 * RcvLoc1 + 4 * RcvLoc2) / 8;
          locy(5) = (3 * RcvLoc1 + 5 * RcvLoc2) / 8;
          locy(6) = (2 * RcvLoc1 + 6 * RcvLoc2) / 8;
          locy(7) = (1 * RcvLoc1 + 7 * RcvLoc2) / 8;
          locy(8) = RcvLoc2;

        } // end of if only recieving one loc data;
        else if (argc - count == 9) {
          double indata[9];
          double BufferData;

          for (int i = 0; i < 9; i++) {
            if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid data " << argv[count]
                     << " for -beamThermal 3D\n";
              return TCL_ERROR;
            }
            indata[i] = BufferData;
            count++;
          }
          locy = Vector(indata, 9);
          // temp1,loc1,temp2,loc2...temp9,loc9
        } // end of if only recieving 9 loc data;

        theLoad = new NodalThermalAction(nodeLoadTag, nodeId, locy, theSeries,
                                         thecrds);
        delete thecrds;
      }
      // end of recieving 9 temp data in external file;
      else {
        opserr << "WARNING NodalThermalAction - invalid dataLen\n";
      }
      // end of definition for different data input length(9 or15)

    }
    // end for detecting source
    else {
      if (argc - count == 4) {
        double t1, t2, locY1, locY2;
        if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid T1 " << argv[count]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 2], &t2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid T1 " << argv[count]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 3], &locY2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }

        theLoad = new NodalThermalAction(nodeLoadTag, nodeId, t1, locY1, t2,
                                         locY2, thecrds);
      }
      // for defining a uniform gradient thermal action
    }
    // end for source or no source
    if (theLoad == 0) {
      opserr << "WARNING NodalLoad - out of memory creating load " << argv[1];
      return TCL_ERROR;
    }
    // get the current pattern tag if no tag given in i/p
    if (userSpecifiedPattern == false) {
      if (theTclLoadPattern == 0) {
        opserr << "WARNING no current load pattern - NodalThermalAction "
               << nodeId;
        return TCL_ERROR;
      } else
        loadPatternTag = theTclLoadPattern->getTag();
    }
#endif
  }
  // end of adding NodalThermalAction -------------end---------Liming,[SIF] 2017

  // start of else block, Liming [SIF]
  else {

    // make sure at least one other argument to contain type of system
    if (argc < (2 + ndf)) {
      opserr << "WARNING bad command - want: load nodeId " << ndf
             << " forces\n";
      return TCL_ERROR;
    }

    // get the id of the node
    int nodeId;
    if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      opserr << "WARNING invalid nodeId: " << argv[1];
      opserr << " - load nodeId " << ndf << " forces\n";
      return TCL_ERROR;
    }

    // get the load vector
    Vector forces(ndf);
    for (int i = 0; i < ndf; i++) {
      double theForce;
      if (Tcl_GetDouble(interp, argv[2 + i], &theForce) != TCL_OK) {
        opserr << "WARNING invalid force " << i + 1 << " - load " << nodeId;
        opserr << " " << ndf << " forces\n";
        return TCL_ERROR;
      } else
        forces(i) = theForce;
    }

    // allow some additional options at end of command
    int endMarker = 2 + ndf;
    while (endMarker != argc) {
      if (strcmp(argv[endMarker], "-const") == 0) {
        // allow user to specify const load
        isLoadConst = true;
      } else if (strcmp(argv[endMarker], "-pattern") == 0) {
        // allow user to specify load pattern other than current
        endMarker++;
        userSpecifiedPattern = true;
        if (endMarker == argc ||
            Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

          opserr << "WARNING invalid patternTag - load " << nodeId << " ";
          opserr << ndf << " forces pattern patterntag\n";
          return TCL_ERROR;
        }
      }
      endMarker++;
    }

    // get the current pattern tag if no tag given in i/p
    if (userSpecifiedPattern == false) {
      if (theTclLoadPattern == 0) {
        opserr << "WARNING no current load pattern - load " << nodeId;
        opserr << " " << ndf << " forces\n";
        return TCL_ERROR;
      } else
        loadPatternTag = theTclLoadPattern->getTag();
    }

    // create the load
    theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
    if (theLoad == 0) {
      opserr << "WARNING ran out of memory for load  - load " << nodeId;
      opserr << " " << ndf << " forces\n";
      return TCL_ERROR;
    }

  } // end of Liming change for nodal thermal action , putting the above block
    // into else{ }

  // add the load to the domain
  if (theTclDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
    opserr << "WARNING TclBasicBuilder - could not add load to domain\n";
    delete theLoad;
    return TCL_ERROR;
  }
  nodeLoadTag++;

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}


int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: mass nodeId " << ndf
           << " mass values\n";
    return TCL_ERROR;
  }

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf, ndf);
  double theMass;
  for (int i = 0; i < ndf; i++) {
    if (Tcl_GetDouble(interp, argv[i + 2], &theMass) != TCL_OK) {
      opserr << "WARNING invalid nodal mass term\n";
      opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
      return TCL_ERROR;
    }
    mass(i, i) = theMass;
  }

  if (theTclDomain->setMass(mass, nodeId) != 0) {
    opserr << "WARNING failed to set mass at node " << nodeId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fix nodeId " << ndf
           << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId - fix nodeId " << ndf
           << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  char buffer[80];
  strcpy(buffer, "");
  // get the fixity condition and add the constraint if fixed
  for (int i = 0; i < ndf; i++) {
    int theFixity;
    if (Tcl_GetInt(interp, argv[2 + i], &theFixity) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i + 1 << " - load " << nodeId;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } else {
      if (theFixity != 0) {

        // create a homogeneous constraint
        SP_Constraint *theSP = new SP_Constraint(nodeId, i, 0.0, true);
        if (theSP == 0) {
          opserr << "WARNING ran out of memory for SP_Constraint ";
          opserr << "fix " << nodeId << " " << ndf << " [0,1] conditions\n";
          return TCL_ERROR;
        }

        // add it to the domain
        if (theTclDomain->addSP_Constraint(theSP) == false) {
          opserr << "WARNING could not add SP_Constraint to domain using fix "
                    "command - node may already be constrained\n";
          sprintf(buffer, "%d ", 0);
          delete theSP;
        } else {
          sprintf(buffer, "%d ", theSP->getTag());
          Tcl_AppendResult(interp, buffer, NULL);
        }
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;
  if (strcmp(argv[argc - 2], "-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixX xLoc " << ndf
           << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the xCrd of nodes to be constrained
  double xLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
    opserr << "WARNING invalid xCrd - fixX xLoc " << ndf
           << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i = 0; i < ndf; i++) {
    if (Tcl_GetInt(interp, argv[2 + i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i + 1 << " - fixX " << xLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2 + ndf], "-tol") == 0)
      if (Tcl_GetDouble(interp, argv[3 + ndf], &tol) != TCL_OK) {
        opserr << "WARNING invalid tol specified - fixX " << xLoc << endln;
        return TCL_ERROR;
      }
  }

  theTclDomain->addSP_Constraint(0, xLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;
  if (strcmp(argv[argc - 2], "-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixY yLoc " << ndf
           << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the yCrd of nodes to be constrained
  double yLoc;
  if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
    opserr << "WARNING invalid yCrd - fixY yLoc " << ndf
           << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i = 0; i < ndf; i++) {
    if (Tcl_GetInt(interp, argv[2 + i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i + 1 << " - fixY " << yLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2 + ndf], "-tol") == 0)
      if (Tcl_GetDouble(interp, argv[3 + ndf], &tol) != TCL_OK) {
        opserr << "WARNING invalid tol specified - fixY " << yLoc << endln;
        return TCL_ERROR;
      }
  }

  theTclDomain->addSP_Constraint(1, yLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;
  if (strcmp(argv[argc - 2], "-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixZ zLoc " << ndf
           << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the yCrd of nodes to be constrained
  double zLoc;
  if (Tcl_GetDouble(interp, argv[1], &zLoc) != TCL_OK) {
    opserr << "WARNING invalid zCrd - fixZ zLoc " << ndf
           << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i = 0; i < ndf; i++) {
    if (Tcl_GetInt(interp, argv[2 + i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i + 1 << " - fixZ " << zLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2 + ndf], "-tol") == 0)
      if (Tcl_GetDouble(interp, argv[3 + ndf], &tol) != TCL_OK) {
        opserr << "WARNING invalid tol specified - fixZ " << zLoc << endln;
        return TCL_ERROR;
      }
  }

  theTclDomain->addSP_Constraint(2, zLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addSP(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";
    return TCL_ERROR;
  }


  // check number of arguments
  if (argc < 4) {
    opserr << "WARNING bad command - want: sp nodeId dofID value";
    return TCL_ERROR;
  }

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId;
  double value;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1]
           << " -  sp nodeId dofID value\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    opserr << "WARNING invalid dofId: " << argv[2] << " -  sp ";
    opserr << nodeId << " dofID value\n";
    return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING invalid value: " << argv[3] << " -  sp ";
    opserr << nodeId << " dofID value\n";
    return TCL_ERROR;
  }

  bool isSpConst = false;
  bool userSpecifiedPattern = false;
  int loadPatternTag = 0; // some pattern that will never be used!

  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker], "-const") == 0) {
      // allow user to specify const load
      isSpConst = true;
    } else if (strcmp(argv[endMarker], "-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      userSpecifiedPattern = true;
      if (endMarker == argc ||
          Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

        opserr << "WARNING invalid patternTag - load " << nodeId << "\n";
        return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (userSpecifiedPattern == false) {
    if (theTclLoadPattern == 0) {
      opserr << "WARNING no current pattern - sp " << nodeId
             << " dofID value\n";
      return TCL_ERROR;
    } else
      loadPatternTag = theTclLoadPattern->getTag();
  }

  LoadPattern *thePattern = theTclDomain->getLoadPattern(loadPatternTag);

  // create a homogeneous constraint
  SP_Constraint *theSP = new SP_Constraint(nodeId, dofId, value, isSpConst);

  if (theSP == 0) {
    opserr << "WARNING ran out of memory for SP_Constraint ";
    opserr << " - sp " << nodeId << " dofID value\n";
    return TCL_ERROR;
  }
  if (theTclDomain->addSP_Constraint(theSP, loadPatternTag) == false) {
    opserr << "WARNING could not add SP_Constraint to domain ";
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addImposedMotionSP(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";
    return TCL_ERROR;
  }


  // check number of arguments
  if (argc < 4) {
    opserr
        << "WARNING bad command - want: imposedMotion nodeId dofID gMotionID\n";
    return TCL_ERROR;
  }

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId, gMotionID;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - imposedMotion nodeId dofID gMotionID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    opserr << "WARNING invalid dofId: " << argv[2] << " -  imposedMotion ";
    opserr << nodeId << " dofID gMotionID\n";
    return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetInt(interp, argv[3], &gMotionID) != TCL_OK) {
    opserr << "WARNING invalid gMotionID: " << argv[3] << " -  imposedMotion ";
    opserr << nodeId << " dofID gMotionID\n";
    return TCL_ERROR;
  }

  bool alt = false;
  if (argc == 5) {
    if (strcmp(argv[4], "-other") == 0)
      alt = true;
  }

  //
  // check valid node & dof
  //

  Node *theNode = theTclDomain->getNode(nodeId);
  if (theNode == 0) {
    opserr << "WARNING invalid node " << argv[2] << " node not found\n ";
    return -1;
  }
  int nDof = theNode->getNumberDOF();
  if (dofId < 0 || dofId >= nDof) {
    opserr << "WARNING invalid dofId: " << argv[2]
           << " dof specified cannot be <= 0 or greater than num dof at nod\n ";
    return -2;
  }

  MultiSupportPattern *thePattern = (MultiSupportPattern *)Tcl_GetAssocData(interp, "theTclMultiSupportPattern", NULL);
  int loadPatternTag = thePattern->getTag();

  // create a new ImposedMotionSP
  SP_Constraint *theSP;
  if (alt == true) {
    theSP = new ImposedMotionSP1(nodeId, dofId, loadPatternTag, gMotionID);
  } else {
    theSP = new ImposedMotionSP(nodeId, dofId, loadPatternTag, gMotionID);
  }

  if (theSP == 0) {
    opserr << "WARNING ran out of memory for ImposedMotionSP ";
    opserr << " -  imposedMotion ";
    opserr << nodeId << " " << dofId++ << " " << gMotionID << endln;
    return TCL_ERROR;
  }
  if (thePattern->addSP_Constraint(theSP) == false) {
    opserr << "WARNING could not add SP_Constraint to pattern ";
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addEqualDOF_MP(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  // Ensure the destructor has not been called
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - equalDOF \n";
    return TCL_ERROR;
  }

  // Check number of arguments
  if (argc < 4) {
    opserr << "WARNING bad command - want: equalDOF RnodeID? CnodeID? DOF1? "
              "DOF2? ...";
    return TCL_ERROR;
  }

  // Read in the node IDs and the DOF
  int RnodeID, CnodeID, dofID;

  if (Tcl_GetInt(interp, argv[1], &RnodeID) != TCL_OK) {
    opserr << "WARNING invalid RnodeID: " << argv[1]
           << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &CnodeID) != TCL_OK) {
    opserr << "WARNING invalid CnodeID: " << argv[2]
           << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
    return TCL_ERROR;
  }

  // The number of DOF to be coupled
  int numDOF = argc - 3;

  // The constraint matrix ... U_c = C_cr * U_r
  Matrix Ccr(numDOF, numDOF);
  Ccr.Zero();

  // The vector containing the retained and constrained DOFs
  ID rcDOF(numDOF);

  int i, j;
  // Read the degrees of freedom which are to be coupled
  for (i = 3, j = 0; i < argc; i++, j++) {
    if (Tcl_GetInt(interp, argv[i], &dofID) != TCL_OK) {
      opserr << "WARNING invalid dofID: " << argv[3]
             << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }

    dofID -= 1; // Decrement for C++ indexing
    if (dofID < 0) {
      opserr << "WARNING invalid dofID: " << argv[i] << " must be >= 1";
      return TCL_ERROR;
    }
    rcDOF(j) = dofID;
    Ccr(j, j) = 1.0;
  }

  // Create the multi-point constraint
  MP_Constraint *theMP = new MP_Constraint(RnodeID, CnodeID, Ccr, rcDOF, rcDOF);
  if (theMP == 0) {
    opserr << "WARNING ran out of memory for equalDOF MP_Constraint ";
    return TCL_ERROR;
  }

  // Add the multi-point constraint to the domain
  if (theTclDomain->addMP_Constraint(theMP) == false) {
    opserr << "WARNING could not add equalDOF MP_Constraint to domain ";
    delete theMP;
    return TCL_ERROR;
  }

  char buffer[80];
  sprintf(buffer, "%d", theMP->getTag());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
TclCommand_addEqualDOF_MP_Mixed(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{
  // Ensure the destructor has not been called
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - equalDOF \n";
    return TCL_ERROR;
  }

  // Check number of arguments
  if (argc < 4) {
    opserr << "WARNING bad command - want: equalDOFmixed RnodeID? CnodeID? "
              "numDOF? RDOF1? CDOF1? ... ...";
    return TCL_ERROR;
  }

  // Read in the node IDs and the DOF
  int RnodeID, CnodeID, dofIDR, dofIDC, numDOF;

  if (Tcl_GetInt(interp, argv[1], &RnodeID) != TCL_OK) {
    opserr << "WARNING invalid RnodeID: " << argv[1]
           << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &CnodeID) != TCL_OK) {
    opserr << "WARNING invalid CnodeID: " << argv[2]
           << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &numDOF) != TCL_OK) {
    opserr << "WARNING invalid numDOF: " << argv[2]
           << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }

  // The number of DOF to be coupled
  //        int numDOF = argc - 3;

  // The constraint matrix ... U_c = C_cr * U_r
  Matrix Ccr(numDOF, numDOF);
  Ccr.Zero();

  // The vector containing the retained and constrained DOFs
  ID rDOF(numDOF);
  ID cDOF(numDOF);

  int i, j, k;
  // Read the degrees of freedom which are to be coupled
  for (i = 4, j = 5, k = 0; k < numDOF; i += 2, j += 2, k++) {
    if (Tcl_GetInt(interp, argv[i], &dofIDR) != TCL_OK) {
      opserr << "WARNING invalid dofID: " << argv[3]
             << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[j], &dofIDC) != TCL_OK) {
      opserr << "WARNING invalid dofID: " << argv[3]
             << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }

    dofIDR -= 1; // Decrement for C++ indexing
    dofIDC -= 1;
    if (dofIDC < 0 || dofIDR < 0) {
      opserr << "WARNING invalid dofID: " << argv[i] << " must be >= 1";
      return TCL_ERROR;
    }
    rDOF(k) = dofIDR;
    cDOF(k) = dofIDC;
    Ccr(k, k) = 1.0;
  }

  // Create the multi-point constraint
  MP_Constraint *theMP = new MP_Constraint(RnodeID, CnodeID, Ccr, cDOF, rDOF);
  if (theMP == 0) {
    opserr << "WARNING ran out of memory for equalDOF MP_Constraint ";
    return TCL_ERROR;
  }

  // Add the multi-point constraint to the domain
  if (theTclDomain->addMP_Constraint(theMP) == false) {
    opserr << "WARNING could not add equalDOF MP_Constraint to domain ";
    delete theMP;
    return TCL_ERROR;
  }

  char buffer[80];
  sprintf(buffer, "%d", theMP->getTag());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

#if 0
int
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  if (argc < 4) {
    opserr << "WARNING rigidLink linkType? rNode? cNode?\n";
    return TCL_ERROR;
  }

  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
    opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read "
              "rNode \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
    opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read "
              "CNode \n";
    return TCL_ERROR;
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1], "-bar") == 0) || (strcmp(argv[1], "bar") == 0)) {
    RigidRod theLink(*theTclDomain, rNode, cNode);
  } else if ((strcmp(argv[1], "-beam") == 0) ||
             (strcmp(argv[1], "beam") == 0)) {
    RigidBeam theLink(*theTclDomain, rNode, cNode);
  } else {
    opserr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link "
              "type (-bar, -beam) \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
#endif


int
TclCommand_addMP(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  opserr << "WARNING - TclCommand_addMP() not yet implemented\n";
  return TCL_OK;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Added by Scott J. Brandenberg, UC Davis, sjbrandenberg@ucdavis.edu
int
TclCommand_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  if (argc < 6 || argc > 7) {
    opserr
        << "WARNING PySimple1Gen file1? file2? file3? file4? file5? <file6?>";
    opserr << "Must have either 5 or 6 arguments." << endln;
  }

  PySimple1Gen *thePySimple1Gen;
  thePySimple1Gen = new PySimple1Gen;

  if (argc == 6)
    thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4],
                                    argv[5]);
  if (argc == 7)
    thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4], argv[5],
                                    argv[6]);

  delete thePySimple1Gen;

  return TCL_OK;
}

int
TclCommand_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  if (argc < 6 || argc > 7) {
    opserr
        << "WARNING TzSimple1Gen file1? file2? file3? file4? file5? <file6?>";
    opserr << "Must have either 5 or 6 arguments." << endln;
  }

  TzSimple1Gen *theTzSimple1Gen;
  theTzSimple1Gen = new TzSimple1Gen;

  if (argc == 6)
    theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4],
                                    argv[5]);
  if (argc == 7)
    theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4], argv[5],
                                    argv[6]);

  delete theTzSimple1Gen;

  return TCL_OK;
}
// End Added by Scott J. Brandenberg
///////////////////////////////////////////////////////////////////////////////////////////////////

// Added by Prishati Raychowdhury (UCSD)
int
TclBasicBuilder_doShallowFoundationGen(ClientData clientData,
                                       Tcl_Interp *interp, int argc,
                                       TCL_Char ** const argv)
{
  if (argc != 5) {
    opserr << "WARNING ShallowFoundationGen FoundationID? ConnectingNode? "
              "InputDataFile? FoundationMatType?";
    opserr << "Must have 4 arguments." << endln;
  }

  ShallowFoundationGen *theShallowFoundationGen;
  theShallowFoundationGen = new ShallowFoundationGen;

  // Checking for error
  int FoundationID;
  int ConnectingNode;
  int FoundationMatType;

  if (Tcl_GetInt(interp, argv[1], &FoundationID) != TCL_OK) {
    opserr << "WARNING invalid FoundationID: " << argv[1]
           << ". ShallowFoundationGen FoundationID? ConnectingNode? "
              "InputDataFile? FoundationMatType? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &ConnectingNode) != TCL_OK) {
    opserr << "WARNING invalid ConnectingNode: " << argv[2]
           << ". ShallowFoundationGen FoundationID? ConnectingNode? "
              "InputDataFile? FoundationMatType? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &FoundationMatType) != TCL_OK) {
    opserr << "WARNING invalid FoundationMatType: " << argv[4]
           << ". ShallowFoundationGen FoundationID? ConnectingNode? "
              "InputDataFile? FoundationMatType? ";
    return TCL_ERROR;
  }

  theShallowFoundationGen->GetShallowFoundation(argv[1], argv[2], argv[3],
                                                argv[4]);
  delete theShallowFoundationGen;

  return TCL_OK;
}
// End PRC

int
TclBasicBuilder_addRemoHFiber(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  return TclCommand_addHFiber(clientData, interp, argc, argv, theTclBuilder);
}

extern int TclBasicBuilderStiffnessDegradationCommand(ClientData clientData,
                                                      Tcl_Interp *interp,
                                                      int argc, TCL_Char ** const argv,
                                                      Domain *theDomain);

int
TclCommand_addStiffnessDegradation(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char ** const argv)
{
  return TclBasicBuilderStiffnessDegradationCommand(clientData, interp, argc,
                                                    argv, theTclDomain);
}

extern int TclBasicBuilderUnloadingRuleCommand(ClientData clientData,
                                               Tcl_Interp *interp, int argc,
                                               TCL_Char ** const argv,
                                               Domain *theDomain);

int
TclCommand_addUnloadingRule(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char ** const argv)
{
  return TclBasicBuilderUnloadingRuleCommand(clientData, interp, argc, argv,
                                             theTclDomain);
}

extern int TclBasicBuilderStrengthDegradationCommand(ClientData clientData,
                                                     Tcl_Interp *interp,
                                                     int argc, TCL_Char ** const argv,
                                                     Domain *theDomain);

int
TclCommand_addStrengthDegradation(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  return TclBasicBuilderStrengthDegradationCommand(clientData, interp, argc,
                                                   argv, theTclDomain);
}

#if 0
/// added by ZHY
extern int TclCommand_UpdateMaterialsCommand(ClientData clientData,
                                             Tcl_Interp *interp, int argc,
                                             TCL_Char ** const argv,
                                             TclBasicBuilder *theTclBuilder,
                                             Domain *theDomain);
int
TclCommand_UpdateMaterials(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  return TclCommand_UpdateMaterialsCommand(clientData, interp, argc, argv,
                                           theTclBuilder, theTclDomain);
}
#endif


extern int TclBasicBuilderFrictionModelCommand(ClientData clienData,
                                               Tcl_Interp *interp, int argc,
                                               TCL_Char ** const argv,
                                               Domain *theDomain);

int
TclCommand_addFrictionModel(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char ** const argv)
{
  return TclBasicBuilderFrictionModelCommand(clientData, interp, argc, argv,
                                             theTclDomain);
}

int
TclCommand_Package(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{

  void *libHandle;
  int (*funcPtr)(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv, Domain *, TclBasicBuilder *);

  const char *funcName = 0;
  int res = -1;

  if (argc == 2) {
    res = getLibraryFunction(argv[1], argv[1], &libHandle, (void **)&funcPtr);
  } else if (argc == 3) {
    res = getLibraryFunction(argv[1], argv[2], &libHandle, (void **)&funcPtr);
  }

  if (res == 0) {
    int result =
        (*funcPtr)(clientData, interp, argc, argv, theTclDomain, theTclBuilder);
  } else {
    opserr << "Error: Could not find function: " << argv[1] << endln;
    return -1;
  }

  return res;
}

