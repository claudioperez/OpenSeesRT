/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.52 $
// $Date: 2010-05-12 20:17:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclBuilder.cpp,v $

// Written: fmk
// Created: 07/99
//
// Description: This file contains the class definition for TclBuilder.
// A TclBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework. currently these elements include:

// What: "@(#) TclBuilder.cpp, revA"

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
#include <TclBuilder.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>

#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>                   //L.Jiang [SIF]
#include <vector>                                    //L.Jiang [SIF]
using std::vector;                                   // L.Jiang [SIF]
#include <SimulationInformation.h>                   //L.Jiang [SIF]
extern SimulationInformation simulationInfo;         // L.Jiang [SIF]
extern const char *getInterpPWD(Tcl_Interp *interp); // L.Jiang [SIF]

#ifdef OPSDEF_ELEMENT_BLOCKND
#include <Block2D.h>
#include <Block3D.h>
#endif // OPSDEF_ELEMENT_BLOCKND
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
#include <CyclicModel.h> //!!
#ifdef OPSDEF_DAMAGE
#include <DamageModel.h> //!!
#endif

#include <FrictionModel.h>

#include <StiffnessDegradation.h>
#include <UnloadingRule.h>
#include <StrengthDegradation.h>
#include <HystereticBackbone.h>
#include <BeamIntegration.h>

////////////////////// gnp adding damping
#include <Element.h>
////////////////////////////////////////////

extern void TCL_OPS_setModelBuilder(TclBuilder *theNewBuilder);
extern int OPS_ResetInput(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char **argv, Domain *domain,
                          TclBuilder *builder);
#include <packages.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Domain *theTclDomain = 0;
static TclBuilder *theTclBuilder = 0;

extern LoadPattern *theTclLoadPattern;
extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
static int eleLoadTag = 0;

//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static int TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);

static int TclCommand_addNode(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char **argv);

static int TclCommand_addElement(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv);

static int TclCommand_mesh(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char **argv);
static int TclCommand_remesh(ClientData clientData, Tcl_Interp *interp,
                             int argc, TCL_Char **argv);
#if defined(OPSDEF_Element_PFEM)
static int TclCommand_backgroundMesh(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char **argv);
#endif // _OPS_Element_PFEM

static int TclCommand_addUniaxialMaterial(ClientData clientData,
                                          Tcl_Interp *interp, int argc,
                                          TCL_Char **argv);

static int TclCommand_addBeamIntegration(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);

static int TclCommand_addLimitCurve(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char **argv);

static int TclCommand_addNDMaterial(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char **argv);

static int TclCommand_addSection(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv);

static int TclCommand_addYieldSurface_BC(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);

static int TclCommand_addYS_EvolutionModel(ClientData clientData,
                                           Tcl_Interp *interp, int argc,
                                           TCL_Char **argv);

static int TclCommand_addYS_PlasticMaterial(ClientData clientData,
                                            Tcl_Interp *interp, int argc,
                                            TCL_Char **argv);

int //!!
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);

#ifdef OPSDEF_DAMAGE
int //!!
TclCommand_addDamageModel(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);
#endif // OPSDEF_DAMAGE

static int TclCommand_addTimeSeries(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char **argv);

static int TclCommand_addPattern(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv);

static int TclCommand_addSeries(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv);

static int TclCommand_addHomogeneousBC(ClientData clientData,
                                       Tcl_Interp *interp, int argc,
                                       TCL_Char **argv);
static int TclCommand_addHomogeneousBC_X(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);
static int TclCommand_addHomogeneousBC_Y(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);
static int TclCommand_addHomogeneousBC_Z(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);
static int TclCommand_addEqualDOF_MP(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char **argv);

static int TclCommand_addEqualDOF_MP_Mixed(ClientData clientData,
                                           Tcl_Interp *interp, int argc,
                                           TCL_Char **argv);

static int TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv);

static int TclCommand_RigidDiaphragm(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char **argv);

static int TclCommand_addMP(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char **argv);

static int TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);

static int TclCommand_addElementalLoad(ClientData clientData,
                                       Tcl_Interp *interp, int argc,
                                       TCL_Char **argv);

static int TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);
static int TclCommand_addSP(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char **argv);

static int TclCommand_addImposedMotionSP(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);
// Added by Scott J. Brandenberg
static int TclCommand_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char **argv);

static int TclCommand_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char **argv);
// End added by SJB

// Added by Prishati Raychowdhury (UCSD)
static int TclBuilder_doShallowFoundationGen(ClientData clientData,
                                                  Tcl_Interp *interp, int argc,
                                                  TCL_Char **argv);
// End PRC

#ifdef OPSDEF_ELEMENT_BLOCKND
static int TclCommand_doBlock2D(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv);

static int TclCommand_doBlock3D(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv);

#endif // OPSDEF_ELEMENT_BLOCKND
static int TclCommand_addRemoPatch(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);

static int TclCommand_addRemoLayer(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);

static int TclCommand_addRemoFiber(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **argv);

// Leo
static int TclBuilder_addRemoHFiber(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);

static int TclCommand_addRemoGeomTransf(ClientData clientData,
                                        Tcl_Interp *interp, int argc,
                                        TCL_Char **argv);

static int TclCommand_addFrictionModel(ClientData clientData,
                                       Tcl_Interp *interp, int argc,
                                       TCL_Char **argv);

static int TclCommand_addStiffnessDegradation(ClientData clientData,
                                              Tcl_Interp *interp, int argc,
                                              TCL_Char **argv);

static int TclCommand_addUnloadingRule(ClientData clientData,
                                       Tcl_Interp *interp, int argc,
                                       TCL_Char **argv);

static int TclCommand_addStrengthDegradation(ClientData clientData,
                                             Tcl_Interp *interp, int argc,
                                             TCL_Char **argv);

static int TclCommand_addHystereticBackbone(ClientData clientData,
                                            Tcl_Interp *interp, int argc,
                                            TCL_Char **argv);

static int TclCommand_addGroundMotion(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char **argv);

/// added by ZHY
static int TclCommand_UpdateMaterialStage(ClientData clientData,
                                          Tcl_Interp *interp, int argc,
                                          TCL_Char **argv);
static int TclCommand_UpdateMaterials(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char **argv);

/// added by ZHY
static int TclCommand_UpdateParameter(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char **argv);

////////////////gnp adding rayleigh //////////////////////////
static int TclCommand_addElementRayleigh(ClientData clientData,
                                         Tcl_Interp *interp, int argc,
                                         TCL_Char **argv);
///////////////////////////////////////////////////////////////

// REMO
extern int TclCommand_addPatch(ClientData clientData, Tcl_Interp *interp,
                               int argc, TCL_Char **argv,
                               TclBuilder *theTclBuilder);

extern int TclCommand_addFiber(ClientData clientData, Tcl_Interp *interp,
                               int argc, TCL_Char **argv,
                               TclBuilder *theTclBuilder);

extern int TclCommand_addHFiber(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv,
                                TclBuilder *theTclBuilder);

extern int TclCommand_addReinfLayer(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char **argv,
                                    TclBuilder *theTclBuilder);

extern int TclCommand_addGeomTransf(ClientData, Tcl_Interp *, int, TCL_Char **,
                                    Domain *, TclBuilder *);

static int TclCommand_Package(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char **argv);


// Added by Alborz Ghofrani - U.Washington
extern int TclCommand_GenerateInterfacePoints(ClientData clientData,
                                              Tcl_Interp *interp, int argc,
                                              TCL_Char **argv);
// End Added by Alborz


//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclBuilder::TclBuilder(Domain &theDomain, int NDM, int NDF) // , Tcl_Interp *interp, int NDM,
                                 // int NDF)
    : ModelBuilder(theDomain), ndm(NDM), ndf(NDF) // , theInterp(interp)
{}

TclBuilder::~TclBuilder(){}

int
TclBuilder::buildFE_Model(void)
{
  return 0;
}

int
TclBuilder::getNDM(void) const
{
  return ndm;
}

int
TclBuilder::getNDF(void) const
{
  return ndf;
}

