// Written: cmp
//
// Description: This file contains the class definition for TclSafeBuilder.
// A TclSafeBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework.
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

#include <UniaxialMaterial.h>
#include <LimitCurve.h>
#include <NDMaterial.h>
// #include <HystereticBackbone.h>
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

// Added by Prishati Raychowdhury  (PRC)
#include <ShallowFoundationGen.h>
//end PRC

#include <YieldSurface_BC.h>
#include <YS_Evolution.h>

#include <HystereticBackbone.h>
#include <BeamIntegration.h>
*/

//////// gnp adding damping ////////////////
#include <Element.h>
////////////////////////////////////////////

/*
extern void TCL_OPS_setModelBuilder(TclSafeBuilder *theNewBuilder);
extern int OPS_ResetInput(ClientData clientData,
                          Tcl_Interp *interp,
                          int cArg,
                          int mArg,
                          TCL_Char **argv,
                          Domain *domain,
                          TclBuilder *builder);
#include <packages.h>
*/


//
// CLASS CONSTRUCTOR & DESTRUCTOR
//
// constructor: the constructor will add certain commands to the interpreter
TclSafeBuilder::TclSafeBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM,
                               int NDF)
    : TclBuilder(theDomain, NDM, NDF), theInterp(interp)
{
  /*
  theSections = new ArrayOfTaggedObjects(32);
  theUniaxialMaterials = new MapOfTaggedObjects();
  theSectionRepresents = new ArrayOfTaggedObjects(32);
  */

  // theYieldSurface_BCs = new ArrayOfTaggedObjects(32);
  // theCycModels = new ArrayOfTaggedObjects(32); //!!
  // theYS_EvolutionModels = new ArrayOfTaggedObjects(32);
  // thePlasticMaterials = new ArrayOfTaggedObjects(32);

  static int ncmd = sizeof(tcl_char_cmds)/sizeof(char_cmd);

  for (int i = 0; i < ncmd; i++)
    Tcl_CreateCommand(interp, 
        tcl_char_cmds[i].name, 
        tcl_char_cmds[i].func, 
        (ClientData)NULL, NULL);
 
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
  // free up memory allocated in the constructor
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
  Tcl_DeleteCommand(theInterp, "node");
  Tcl_DeleteCommand(theInterp, "element");
  Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  Tcl_DeleteCommand(theInterp, "nDMaterial");
  Tcl_DeleteCommand(theInterp, "section");
  Tcl_DeleteCommand(theInterp, "pattern");
  Tcl_DeleteCommand(theInterp, "timeSeries");
  Tcl_DeleteCommand(theInterp, "load");
}

//
// CLASS METHODS
/*
int TclSafeBuilder::getNDM(void) const {return ndm;}
int TclSafeBuilder::getNDF(void) const {return ndf;}
int TclSafeBuilder::buildFE_Model(void) {return 0;}
*/
int TclSafeBuilder::incrNodalLoadTag(void){return ++nodeLoadTag;};
int TclSafeBuilder::decrNodalLoadTag(void){return --nodeLoadTag;};
int TclSafeBuilder::getNodalLoadTag(void) const {return   nodeLoadTag;};

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

int
TclSafeBuilder::addUniaxialMaterial(UniaxialMaterial *mat)
{
  return this->addUniaxialMaterial(*mat);
}

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
  G3_Runtime *rt = G3_getRuntime(interp);
  TclSafeBuilder *theTclBuilder = G3_getSafeBuilder(rt);
  Domain *theTclDomain = G3_getDomain(rt);

  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }


  int ndm = G3_getNDM(rt);
  int ndf = G3_getNDF(rt);

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
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid YCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid YCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      opserr << "WARNING invalid ZCoordinate\n";
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

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

int
TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  TclSafeBuilder *theTclBuilder = G3_getSafeBuilder(rt);
  Domain *theTclDomain = G3_getDomain(rt);
  int nodeLoadTag = theTclBuilder->getNodalLoadTag();
  LoadPattern *theTclLoadPattern = theTclBuilder->getEnclosingPattern();
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";
    return TCL_ERROR;
  }

  int ndf = argc - 2;
  NodalLoad *theLoad = 0;

  bool isLoadConst = false;
  bool userSpecifiedPattern = false;
  int loadPatternTag = 0;
  // The above definition are moved forward for the use in both cases

  //-------------Adding Proc for NodalThermalAction, By Liming Jiang, [SIF] 2017
  if ((strcmp(argv[2], "-NodalThermal") == 0) ||
      (strcmp(argv[2], "-nodalThermal") == 0)) {

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
      // simulationInfo.addInputFile(argv[count], pwd);
      TimeSeries *theSeries;
      // default num of temperature input for nodal ThermalAction;
      int dataLen = 9; 

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
      // end of definition for different data input length(9or15)
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
        opserr << "WARNING no current load pattern -NodalThermalAction "
               << nodeId;
        return TCL_ERROR;
      } else
        loadPatternTag = theTclLoadPattern->getTag();
    }
  }
  // end of adding NodalThermalAction -------------end--------- Liming,[SIF]2017

      // start of else block, Liming [SIF]
      else
  {

    // make sure at least one other argument to contain type of system
    if (argc < (2 + ndf)) {
      opserr << "WARNING bad command - want: load nodeId " << ndf << "forces\n";
      printCommand(argc, argv);
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
        opserr << "WARNING invalid force " << i + 1 << " - load" << nodeId;
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
    opserr << "WARNING TclSafeBuilder - could not add load to domain\n";
    printCommand(argc, argv);
    delete theLoad;
    return TCL_ERROR;
  }
  theTclBuilder->incrNodalLoadTag();

  // if get here we have sucessfully created the load and added it to the domain
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
TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp, int argc,
TCL_Char **argv)

{
  return TclSafeBuilderParameterCommand(clientData, interp,
                                         argc, argv, theTclDomain,
theTclBuilder);
}

*/
class TclBasicBuilder;
extern int TclBasicBuilderElementCommand(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv, Domain *theDomain,
                                         TclBasicBuilder *theTclBuilder);

int
TclCommand_addElement(ClientData clientData, Tcl_Interp *interp, int argc,
                      TCL_Char **argv)

{
  G3_Runtime *rt = G3_getRuntime(interp);
  TclSafeBuilder *theTclBuilder = G3_getSafeBuilder(rt);
  Domain *theTclDomain = G3_getDomain(rt);
  return TclBasicBuilderElementCommand(clientData, interp, argc, argv,
                                       theTclDomain,
                                       (TclBasicBuilder *)theTclBuilder);
}

extern int TclBasicBuilderUniaxialMaterialCommand(ClientData, Tcl_Interp *, int,
                                                  TCL_Char **, Domain *);

extern int TclSafeBuilderUniaxialCommand(ClientData, Tcl_Interp *, int,
                                         TCL_Char **, Domain *);
int
TclCommand_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
                               int argc, TCL_Char **argv)
{
  int stat = TclSafeBuilderUniaxialCommand(clientData, interp, argc, argv,
                                       G3_getDomain(G3_getRuntime(interp)));
  return stat;
}

/*
extern int
TclSafeBuilderNDMaterialCommand (ClientData clienData, Tcl_Interp
*interp, int argc, TCL_Char **argv, TclSafeBuilder *theTclBuilder);

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

  TimeSeries *theSeries = TclTimeSeriesCommand(clientData, interp, argc - 1,
                                               &argv[1], theTclDomain);

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
    G3_Runtime *rt = G3_getRuntime(interp);
    TclBuilder *theTclBuilder = G3_getModelBuilder(rt);
    Domain     *theTclDomain = G3_getDomain(rt);

  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
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



// Added by Prishati Raychowdhury (UCSD)
int
TclSafeBuilder_doShallowFoundationGen(ClientData clientData, Tcl_Interp
*interp, int argc, TCL_Char **argv)
{
        if(argc != 5){
                opserr << "WARNING ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType?"; opserr << "Must have 4 arguments." << endln;
        }

        ShallowFoundationGen *theShallowFoundationGen;
        theShallowFoundationGen = new ShallowFoundationGen;


      // Checking for error
        int FoundationID; int ConnectingNode; int FoundationMatType;

        if (Tcl_GetInt(interp, argv[1], &FoundationID) != TCL_OK) {
          opserr << "WARNING invalid FoundationID: " << argv[1]
               << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? ";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2], &ConnectingNode) != TCL_OK) {
          opserr << "WARNING invalid ConnectingNode: " << argv[2]
               << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? ";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[4], &FoundationMatType) != TCL_OK) {
          opserr << "WARNING invalid FoundationMatType: " << argv[4]
               << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? "; return TCL_ERROR;
        }

        theShallowFoundationGen->GetShallowFoundation(argv[1], argv[2], argv[3], argv[4]); delete theShallowFoundationGen;

        return TCL_OK;
}
// End PRC
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

// Add a new UniaxialMaterial to the model runtime
int
TclSafeBuilder::addUniaxialMaterial(const std::string &name, UniaxialMaterial &instance)
{
  m_UniaxialMaterialMap[name] = &instance;
  return 1;
}

// Add a new UniaxialMaterial to the model runtime
int
TclSafeBuilder::addUniaxialMaterial(UniaxialMaterial &instance)
{
  const std::string &name = std::to_string(instance.getTag());
  m_UniaxialMaterialMap[name] = &instance;
/*
  opserr << "WARNING (ModelBuilder) Failed to add UniaxialMaterial \n"
         << "         with tag '" << name.c_str() << "' to model.\n";
*/
  return 1;
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

