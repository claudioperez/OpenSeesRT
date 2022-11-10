/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */


// Written: cmp
//
// Description: This file contains the class definition for TclSafeBuilder.
// A TclSafeBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework.
#include <assert.h>
#include <modeling/commands.h>
#include <g3_api.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>

#include <RigidRod.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>

#include <CrdTransf.h>

#include <SectionForceDeformation.h>
#include <SectionRepres.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <TclSafeBuilder.h>
#include <MultiSupportPattern.h>

#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h> //L.Jiang [SIF]

/*
#include <SimulationInformation.h>				//L.Jiang [SIF]
extern SimulationInformation simulationInfo;		//L.Jiang [SIF]
*/
extern const char * getInterpPWD(Tcl_Interp *interp);  //L.Jiang [SIF]

/*--------------------------------------------------------------------

#include <YieldSurface_BC.h>
#include <YS_Evolution.h>

#include <HystereticBackbone.h>
#include <BeamIntegration.h>
*/

#include <Element.h>

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//
// constructor: the constructor will add certain commands to the interpreter
TclSafeBuilder::TclSafeBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM,
                               int NDF)
    : TclBuilder(theDomain, NDM, NDF), theInterp(interp)
{
  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);

  for (int i = 0; i < ncmd; i++)
    Tcl_CreateCommand(interp, 
        tcl_char_cmds[i].name, 
        tcl_char_cmds[i].func, 
        (ClientData) this, nullptr);
 
  theTclBuilder = this;
  theTclDomain = &theDomain;
  tclEnclosingPattern = 0;
  // theTclMultiSupportPattern = 0;

  nodeLoadTag = 0;
  eleArgStart = 0;
  m_runtime = G3_getRuntime(interp);

  Tcl_SetAssocData(interp, "OPS::theTclBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theTclSafeBuilder", NULL, (ClientData)this);
  G3_setDomain(m_runtime, &theDomain);
  Tcl_SetAssocData(interp, "OPS::theTclDomain", NULL, (ClientData)&theDomain);
}

TclSafeBuilder::~TclSafeBuilder()
{
  /*
  OPS_clearAllTimeSeries();
  OPS_clearAllUniaxialMaterial();
  OPS_clearAllNDMaterial();
  OPS_clearAllSectionForceDeformation();
  OPS_clearAllCrdTransf();
  */
  // OPS_clearAllHystereticBackbone();
  // OPS_clearAllFrictionModel();
  // OPS_clearAllLimitCurve();
  // OPS_clearAllDamageModel();
  // theYieldSurface_BCs->clearAll();
  // theYS_EvolutionModels->clearAll();
  // thePlasticMaterials->clearAll();
  // theCycModels->clearAll(); //!!

/*
  theSections->clearAll();
  theSectionRepresents->clearAll();
  delete theSections;
  delete theSectionRepresents;
*/

  // set the pointers to 0
  theTclDomain = 0;
  theTclBuilder = 0;
  tclEnclosingPattern = 0;
  // theTclMultiSupportPattern = 0;
  /* TCL_OPS_setModelBuilder(0); */

  // may possibly invoke Tcl_DeleteCommand() later
  // Tcl_DeleteCommand(theInterp, "node");
  // Tcl_DeleteCommand(theInterp, "element");
  // Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  // Tcl_DeleteCommand(theInterp, "nDMaterial");
  // Tcl_DeleteCommand(theInterp, "section");
  // Tcl_DeleteCommand(theInterp, "pattern");
  // Tcl_DeleteCommand(theInterp, "timeSeries");
  // Tcl_DeleteCommand(theInterp, "load");
}

//
// CLASS METHODS
//
void
TclSafeBuilder::letClobber(bool let_clobber) {no_clobber = !let_clobber;};

bool
TclSafeBuilder::canClobber() {return !no_clobber;};

int TclSafeBuilder::incrNodalLoadTag(void){return ++nodeLoadTag;};
int TclSafeBuilder::decrNodalLoadTag(void){return --nodeLoadTag;};
int TclSafeBuilder::getNodalLoadTag(void) {return   nodeLoadTag;};

int
TclSafeBuilder::addSP_Constraint(int axisDirn, double axisValue, const ID &fixityCodes, double tol)
{
  return theTclDomain->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
}

LoadPattern *
TclSafeBuilder::getEnclosingPattern(void) const {return tclEnclosingPattern;};

int
TclSafeBuilder::setEnclosingPattern(LoadPattern* pat){
  tclEnclosingPattern = pat;
  return 1;
};

Domain *
TclSafeBuilder::getDomain(void) const {return theTclDomain;}

TclSafeBuilder *
TclSafeBuilder::getBuilder(void) const {return theTclBuilder;}

TimeSeries *
TclSafeBuilder::getTimeSeries(const std::string &name)
{
  TimeSeries *series = m_TimeSeriesMap.at(name);
  if (series)
    return series->getCopy();
  else
    return 0;
}

int
TclSafeBuilder::addTimeSeries(const std::string &name, TimeSeries *series)
{
  m_TimeSeriesMap[name] = series;
  return 1;
}

int
TclSafeBuilder::addTimeSeries(TimeSeries *series)
{
  const std::string &name = std::to_string(series->getTag());
  m_TimeSeriesMap[name] = series;
  return 1;
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//
static void
printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i = 0; i < argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
}

static int
TclCommand_addNode(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv)
{
  assert(clientData != nullptr);

  TclSafeBuilder *theTclBuilder = (TclSafeBuilder*)clientData;

  Domain *theTclDomain = theTclBuilder->getDomain();

  // TclSafeBuilder *builder = (TclSafeBuilder*)clientData;

  // ensure the destructor has not been called -
  if (theTclBuilder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2 + ndm) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;

  // get the nodal id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeTag\n";
    opserr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) {
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  }

  else if (ndm == 2) {
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid 1st coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid 2nd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid 1st coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid 2nd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      opserr << "WARNING invalid 3rd coordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  } else {
    opserr << "WARNING invalid ndm\n";
    opserr << "node: " << nodeId << endln;
    ;
    return TCL_ERROR;
  }

  // check for -ndf override option
  int currentArg = 2 + ndm;
  if (currentArg < argc && strcmp(argv[currentArg], "-ndf") == 0) {
    if (Tcl_GetInt(interp, argv[currentArg + 1], &ndf) != TCL_OK) {
      opserr << "WARNING invalid nodal ndf given for node " << nodeId << endln;
      return TCL_ERROR;
    }
    currentArg += 2;
  }

  //
  // create the node
  //
  if (ndm == 1)
    theNode = new Node(nodeId, ndf, xLoc);
  else if (ndm == 2)
    theNode = new Node(nodeId, ndf, xLoc, yLoc);
  else
    theNode = new Node(nodeId, ndf, xLoc, yLoc, zLoc);

  //
  // add the node to the domain
  //
  if (theTclDomain->addNode(theNode) == false) {
    opserr << "WARNING failed to add node to the domain\n";
    opserr << "node: " << nodeId << endln;
    delete theNode; // otherwise memory leak
    return TCL_ERROR;
  }

  while (currentArg < argc) {
    if (strcmp(argv[currentArg], "-mass") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal mass terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Matrix mass(ndf, ndf);
      double theMass;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theMass) != TCL_OK) {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        mass(i, i) = theMass;
      }
      theNode->setMass(mass);
    } else if (strcmp(argv[currentArg], "-dispLoc") == 0) {
      currentArg++;
      if (argc < currentArg + ndm) {
        opserr << "WARNING incorrect number of nodal display location terms, "
                  "need ndm\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector displayLoc(ndm);
      double theCrd;
      for (int i = 0; i < ndm; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theCrd) != TCL_OK) {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        displayLoc(i) = theCrd;
      }
      theNode->setDisplayCrds(displayLoc);

    } else if (strcmp(argv[currentArg], "-disp") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal disp terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << "WARNING invalid nodal disp term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialDisp(disp);
      theNode->commitState();

    } else if (strcmp(argv[currentArg], "-vel") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << "WARNING incorrect number of nodal vel terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; i++) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << "WARNING invalid nodal vel term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << endln;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialVel(disp);
      theNode->commitState();

    } else
      currentArg++;
  }

  // if we get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

/*
extern int
TclSafeBuilderParameterCommand(ClientData clientData,
                                Tcl_Interp *interp, int argc,
                                TCL_Char **argv,
                                Domain *theDomain,
                                TclSafeBuilder *theTclBuilder);

int
TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TclSafeBuilderParameterCommand(clientData, interp,
                                         argc, argv, theTclDomain, theTclBuilder);
}

*/


/*
extern int
TclSafeBuilderNDMaterialCommand(ClientData clienData, Tcl_Interp *interp, int argc, 
                                TCL_Char **argv, TclSafeBuilder *theTclBuilder);

int
TclCommand_addNDMaterial(ClientData clientData, Tcl_Interp *interp,
                            int argc,    TCL_Char **argv)

{
  return TclSafeBuilderNDMaterialCommand(clientData, interp,
                                          argc, argv, theTclBuilder);
}
*/

extern int TclPatternCommand(ClientData clientData, Tcl_Interp *interp,
                             int argc, TCL_Char **argv, Domain *theDomain);

static int
TclCommand_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,
                      TCL_Char **argv)
{
  TclSafeBuilder *theTclBuilder =
      (TclSafeBuilder *)Tcl_GetAssocData(interp, "OPS::theTclBuilder", NULL);
  Domain *theTclDomain = theTclBuilder->getDomain();
  return TclPatternCommand(clientData, interp, argc, argv, theTclDomain);
}

extern TimeSeries *TclTimeSeriesCommand(ClientData clientData,
                                        Tcl_Interp *interp, int argc,
                                        TCL_Char **argv, Domain *theDomain);

static int
TclCommand_addTimeSeries(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char **argv)
{
  TclSafeBuilder *theTclBuilder = (TclSafeBuilder *)Tcl_GetAssocData(
      interp,"OPS::theTclSafeBuilder", NULL);
  Domain *theTclDomain = theTclBuilder->getDomain();

  TimeSeries *theSeries = TclTimeSeriesCommand(clientData, interp, argc - 1, &argv[1], theTclDomain);

  if (theSeries != 0) {
    if (theTclBuilder->addTimeSeries(argv[2], theSeries))
      return TCL_OK;
    else
      return TCL_ERROR;
  }
  return TCL_ERROR;
}

int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char **argv)
{
  assert(clientData != nullptr);

  TclSafeBuilder *theTclBuilder = (TclSafeBuilder*)clientData;

  Domain *theTclDomain = theTclBuilder->getDomain();

  if (theTclBuilder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed - load \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: mass nodeId " << ndf << " mass values\n"; 
    printCommand(argc, argv); 
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
  Matrix mass(ndf,ndf);
  double theMass;
  for (int i=0; i<ndf; i++)
  {
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK)
     {
          opserr << "WARNING invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
          return TCL_ERROR;
      }
      mass(i,i) = theMass;
  }

  if (theTclDomain->setMass(mass, nodeId) != 0) {
    opserr << "WARNING failed to set mass at node " << nodeId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}




/*

int
TclCommand_addMP(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char **argv)
{
  opserr << "WARNING - TclCommand_addMP() not yet implemented\n";
  return TCL_OK;
}

*/



/*
int
TclSafeBuilder_addRemoHFiber(ClientData clientData, Tcl_Interp *interp,
int argc, TCL_Char **argv)
{
  return TclCommand_addHFiber(clientData, interp, argc,argv,theTclBuilder);
}


/// added by ZHY
extern int
TclSafeBuilderUpdateMaterialStageCommand(ClientData clientData,
                                          Tcl_Interp *interp,
                                          int argc,
                                          TCL_Char **argv,
                                          TclSafeBuilder *theTclBuilder,
                                          Domain *theDomain);
int
TclCommand_UpdateMaterialStage(ClientData clientData,
                                    Tcl_Interp *interp,
                                    int argc,
                                    TCL_Char **argv)
{
  return TclSafeBuilderUpdateMaterialStageCommand(clientData, interp,
                                                   argc, argv, theTclBuilder,
theTclDomain);
}

/// added by ZHY
extern int
TclCommand_UpdateMaterialsCommand(ClientData clientData,
                                  Tcl_Interp *interp,
                                  int argc,
                                  TCL_Char **argv,
                                  TclSafeBuilder *theTclBuilder,
                                  Domain *theDomain);
static int
TclCommand_UpdateMaterials(ClientData clientData,
                           Tcl_Interp *interp,
                           int argc,
                           TCL_Char **argv)
{
  TclSafeBuilder *theTclBuilder =
      (TclSafeBuilder *)Tcl_GetAssocData(interp, "OPS::theTclBuilder", NULL);
  return TclCommand_UpdateMaterialsCommand(clientData, interp,
                                           argc, argv, theTclBuilder, theTclDomain);
}

/// added by ZHY
extern int
TclSafeBuilderUpdateParameterCommand(ClientData clientData,
                                          Tcl_Interp *interp,
                                          int argc,
                                          TCL_Char **argv,
                                          TclSafeBuilder *theTclBuilder); 
int TclCommand_UpdateParameter(ClientData clientData,
                                    Tcl_Interp *interp,
                                    int argc,
                                    TCL_Char **argv)
{
  return TclSafeBuilderUpdateParameterCommand(clientData, interp,
                                       argc, argv, theTclBuilder);
}


*/


//
// BEGIN AUTGEN
//

//
// SectionForceDeformation Operations
//

// Retrieve a SectionForceDeformation instance from the model
// runtime
SectionForceDeformation*
TclSafeBuilder::getSection(const std::string &name)
{
  SectionForceDeformation *instance = m_SectionForceDeformationMap.at(name);
  if (instance) {
    return instance->getCopy();
  } else {
    return nullptr;
  }
}

SectionForceDeformation*
TclSafeBuilder::getSection(int tag)
{
  const std::string &name = std::to_string(tag);
  return this->getSection(name);
}

// Add a new SectionForceDeformation to the model runtime
int
TclSafeBuilder::addSection(const std::string &name, SectionForceDeformation &instance)
{
  m_SectionForceDeformationMap[name] = &instance;
  return 1;
}

// Add a new SectionForceDeformation to the model runtime
int
TclSafeBuilder::addSection(SectionForceDeformation &instance)
{
  const std::string &name = std::to_string(instance.getTag());
  m_SectionForceDeformationMap[name] = &instance;
/*
  opserr << "WARNING (ModelBuilder) Failed to add SectionForceDeformation \n"
         << "         with tag '" << name.c_str() << "' to model.\n";
*/
  return 1;
}

//
// SectionRepres Operations
//

// Retrieve a SectionRepres instance from the model
// runtime
SectionRepres*
TclSafeBuilder::getSectionRepres(const std::string &name)
{
  SectionRepres *instance = m_SectionRepresMap.at(name);
  if (instance) {
    return instance;
  } else {
    return nullptr;
  }
}

SectionRepres*
TclSafeBuilder::getSectionRepres(int tag)
{
  const std::string &name = std::to_string(tag);
  return this->getSectionRepres(name);
}

// Add a new SectionRepres to the model runtime
int
TclSafeBuilder::addSectionRepres(const std::string &name, SectionRepres &instance)
{
  m_SectionRepresMap[name] = &instance;
  return 1;
}

// Add a new SectionRepres to the model runtime
int
TclSafeBuilder::addSectionRepres(SectionRepres &instance)
{
  const std::string &name = std::to_string(instance.getTag());
  m_SectionRepresMap[name] = &instance;
/*
  opserr << "WARNING (ModelBuilder) Failed to add SectionRepres \n"
         << "         with tag '" << name.c_str() << "' to model.\n";
*/
  return 1;
}

//
// NDMaterial Operations
//

// Retrieve a NDMaterial instance from the model
// runtime
NDMaterial*
TclSafeBuilder::getNDMaterial(const std::string &name)
{
  NDMaterial *instance = m_NDMaterialMap.at(name);
  if (instance) {
    return instance->getCopy();
  } else {
    return nullptr;
  }
}

NDMaterial*
TclSafeBuilder::getNDMaterial(int tag)
{
  const std::string &name = std::to_string(tag);
  return this->getNDMaterial(name);
}

// Add a new NDMaterial to the model runtime
int
TclSafeBuilder::addNDMaterial(const std::string &name, NDMaterial &instance)
{
  m_NDMaterialMap[name] = &instance;
  return 1;
}

// Add a new NDMaterial to the model runtime
int
TclSafeBuilder::addNDMaterial(NDMaterial &instance)
{
  const std::string &name = std::to_string(instance.getTag());
  m_NDMaterialMap[name] = &instance;
/*
  opserr << "WARNING (ModelBuilder) Failed to add NDMaterial \n"
         << "         with tag '" << name.c_str() << "' to model.\n";
*/
  return 1;
}

//
// UniaxialMaterial Operations
//

// Retrieve a UniaxialMaterial instance from the model
// runtime
UniaxialMaterial*
TclSafeBuilder::getUniaxialMaterial(const std::string &name)
{
  UniaxialMaterial *instance = m_UniaxialMaterialMap.at(name);
  if (instance) {
    return instance->getCopy();
  } else {
    return nullptr;
  }
}

UniaxialMaterial*
TclSafeBuilder::getUniaxialMaterial(int tag)
{
  const std::string &name = std::to_string(tag);
  return this->getUniaxialMaterial(name);
}

int
TclSafeBuilder::addUniaxialMaterial(UniaxialMaterial *mat)
{
  return this->addUniaxialMaterial(*mat);
}

// Add a new UniaxialMaterial to the model runtime
int
TclSafeBuilder::addUniaxialMaterial(UniaxialMaterial &instance)
{
  const std::string &name = std::to_string(instance.getTag());
  return this->addUniaxialMaterial(name, instance);
  // m_UniaxialMaterialMap[name] = &instance;
  // return 1;
}

// Add a new UniaxialMaterial to the model runtime
int
TclSafeBuilder::addUniaxialMaterial(const std::string &name, UniaxialMaterial &instance)
{
  if (!canClobber() && (m_UniaxialMaterialMap.find(name) != m_UniaxialMaterialMap.end())) {
    return -1;
  }
  m_UniaxialMaterialMap[name] = &instance;
  return TCL_OK;
}


HystereticBackbone*
TclSafeBuilder::getHystereticBackbone(const std::string &name)
{
  HystereticBackbone *instance = m_HystereticBackboneMap.at(name);
  if (instance) {
    return instance;
  } else {
    return nullptr;
  }
}

// Add a new HystereticBackbone to the model runtime
int
TclSafeBuilder::addHystereticBackbone(const std::string &name, HystereticBackbone &instance)
{
  m_HystereticBackboneMap[name] = &instance;
  return 1;
}



//
// CrdTransf Operations
//

// Retrieve a CrdTransf instance from the model
// runtime
CrdTransf*
TclSafeBuilder::getCrdTransf(const std::string &name)
{
  CrdTransf *instance = m_CrdTransfMap.at(name);
  if (instance) {
    return instance;
  } else {
    return nullptr;
  }
}

CrdTransf*
TclSafeBuilder::getCrdTransf(int tag)
{
  const std::string &name = std::to_string(tag);
  return this->getCrdTransf(name);
}

// Add a new CrdTransf to the model runtime
int
TclSafeBuilder::addCrdTransf(const std::string name, CrdTransf *instance)
{
  // m_CrdTransfMap[name] = instance;
  m_CrdTransfMap.insert({name, instance});
  return 1;
}

// Add a new CrdTransf to the model runtime
int
TclSafeBuilder::addCrdTransf(CrdTransf *instance)
{
  const key_t name = std::to_string(instance->getTag());
  // m_CrdTransfMap[name]std::stringnstance;
  // m_CrdTransfMap.insert(std::make_pair<key_t,CrdTransf*>(std::move(name), instance);
  return this->addCrdTransf(name, instance);
}

