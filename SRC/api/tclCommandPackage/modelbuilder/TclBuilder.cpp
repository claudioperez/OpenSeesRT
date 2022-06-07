

// written: cmp
//
// Description: This file contains the class definition for TclBuilder.


#include <stdlib.h>
#include <string.h>

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

extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
static int eleLoadTag = 0;


TclBuilder::TclBuilder(Domain &theDomain, int NDM, int NDF)
    : ModelBuilder(theDomain), ndm(NDM), ndf(NDF)
{}

TclBuilder::~TclBuilder(){}

int
TclBuilder::buildFE_Model(void) {return 0;}

int
TclBuilder::getNDM(void) const {return ndm;}

int
TclBuilder::getNDF(void) const {return ndf;}

LoadPattern*
TclBuilder::getCurrentLoadPattern(void) {return m_current_load_pattern;}

