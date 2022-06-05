//
// cmp
//
#include <tcl.h>

typedef int (TclCharFn)(ClientData, Tcl_Interp*, int, const char**);
// typedef int (TclObjFn)(ClientData,  Tcl_Interp*, int, Tcl_Obj**);

static TclCharFn  TclCommand_addNode;
static TclCharFn  TclCommand_addPattern;
static TclCharFn  TclCommand_addTimeSeries;
static TclCharFn  TclCommand_addNodalMass;
extern TclCharFn  TclCommand_addGeomTransf;

extern TclCharFn  TclCommand_addElement;
extern TclCharFn  TclCommand_doBlock2D;
extern TclCharFn  TclCommand_doBlock3D;
extern TclCharFn  TclCommand_addUniaxialMaterial;
extern TclCharFn  TclCommand_addHystereticBackbone;
extern TclCharFn  TclCommand_addSection;
extern TclCharFn  TclCommand_addPatch;
extern TclCharFn  TclCommand_addReinfLayer;
// extern TclCharFn  TclCommand_addRemoFiber;
extern TclCharFn  TclCommand_addFiber;
// Constraints
extern TclCharFn  TclCommand_addHomogeneousBC;
extern TclCharFn  TclCommand_addEqualDOF_MP;
// Loads
extern TclCharFn  TclCommand_addNodalLoad;
TclCharFn TclCommand_addElementalLoad;

// Damping
TclCharFn modalDamping;
TclCharFn modalDampingQ;

TclCharFn TclCommand_addParameter;
TclCharFn TclCommand_mesh;
TclCharFn TclCommand_remesh;
TclCharFn TclCommand_backgroundMesh; 
TclCharFn TclCommand_addBeamIntegration;

TclCharFn TclCommand_addLimitCurve;
TclCharFn TclCommand_addNDMaterial;
TclCharFn TclCommand_addSeries;

// Constraints
TclCharFn TclCommand_addHomogeneousBC_X;
TclCharFn TclCommand_addHomogeneousBC_Y; 
TclCharFn TclCommand_addHomogeneousBC_Z;
TclCharFn TclCommand_addEqualDOF_MP_Mixed;
TclCharFn TclCommand_addMP;
TclCharFn TclCommand_addSP;
TclCharFn TclCommand_RigidLink;
TclCharFn TclCommand_addImposedMotionSP;
TclCharFn TclCommand_addGroundMotion;
TclCharFn TclCommand_RigidDiaphragm;


struct char_cmd {
  const char* name;
  TclCharFn*  func;
  bool was_added = false;

}  const tcl_char_cmds[] =  {

  {"node",             TclCommand_addNode},
  {"mass",             TclCommand_addNodalMass},
  {"element",          TclCommand_addElement},

// Materials & sections
  {"uniaxialMaterial", TclCommand_addUniaxialMaterial},
  {"nDMaterial",       TclCommand_addNDMaterial},

  {"section",          TclCommand_addSection},
  {"patch",            TclCommand_addPatch},
  {"fiber",            TclCommand_addFiber},
  {"layer",            TclCommand_addReinfLayer},

  {"geomTransf",       TclCommand_addGeomTransf},

  {"load",             TclCommand_addNodalLoad},
  {"pattern",          TclCommand_addPattern},
  {"timeSeries",       TclCommand_addTimeSeries},

  {"fix",                  TclCommand_addHomogeneousBC},
  {"fixX",                 TclCommand_addHomogeneousBC_X},
  {"fixY",                 TclCommand_addHomogeneousBC_Y},
  {"fixZ",                 TclCommand_addHomogeneousBC_Z},
  {"equalDOF",             TclCommand_addEqualDOF_MP},
  {"rigidLink",            &TclCommand_RigidLink},

  {"sp",                   TclCommand_addSP},
  {"groundMotion",         TclCommand_addGroundMotion},
  {"imposedMotion",        TclCommand_addImposedMotionSP},
  {"imposedSupportMotion", TclCommand_addImposedMotionSP},


  {"modalDamping",     modalDamping},
  {"modalDampingQ",    modalDampingQ},


/*
  {"beamIntegration",  TclCommand_addBeamIntegration},
  {"eleLoad",          TclCommand_addElementalLoad},
*/

/*
  {"mp",                   TclCommand_addMP},

  {"block2D",              TclCommand_doBlock2D},
  {"block3D",              TclCommand_doBlock3D},

  {"equalDOF_Mixed",       TclCommand_addEqualDOF_MP_Mixed},
  {"rigidDiaphragm",       &TclCommand_RigidDiaphragm},
  {"PySimple1Gen",         TclCommand_doPySimple1Gen},
  {"TzSimple1Gen",         TclCommand_doTzSimple1Gen},
  {"ShallowFoundationGen", TclSafeBuilder_doShallowFoundationGen},
  {"Hfiber",               TclSafeBuilder_addRemoHFiber},

#if defined(OPSDEF_Element_PFEM)
  {"mesh",             TclCommand_mesh},
  {"remesh",           TclCommand_remesh},
  {"background",      &TclCommand_backgroundMesh},
#endif // OPSDEF_Element_PFEM
*/

// OTHER OBJECT TYPES
  {"hystereticBackbone",   TclCommand_addHystereticBackbone},
  {          "backbone",   TclCommand_addHystereticBackbone},

/*
  {"yieldSurface_BC",      TclCommand_addYieldSurface_BC},
  {"ysEvolutionModel",     TclCommand_addYS_EvolutionModel},
  {"plasticMaterial",      TclCommand_addYS_PlasticMaterial},
  {"cyclicModel",          TclCommand_addCyclicModel},
  {"limitCurve",           TclCommand_addLimitCurve},
  {"damageModel",          TclCommand_addDamageModel},
  {"frictionModel",        TclCommand_addFrictionModel},
  {"stiffnessDegradation", TclCommand_addStiffnessDegradation},
  {"unloadingRule",        TclCommand_addUnloadingRule},
  {"strengthDegradation",  TclCommand_addStrengthDegradation},
  {"loadPackage",          TclCommand_Package},
*/

/*
  {"parameter",            TclCommand_addParameter},
  {"addToParameter",       TclCommand_addParameter},
  {"updateParameter",      TclCommand_addParameter},
*/

/*
// new command for elast2plast in Multi-yield
// plasticity, by ZHY
  {"updateMaterialStage", TclCommand_UpdateMaterialStage},
  {"updateMaterials",     TclCommand_UpdateMaterials},
*/

/*
// new command for updating properties of
// soil materials, by ZHY
   {"updateParameter", TclCommand_UpdateParameter},
*/

};



TclCharFn TclCommand_Package;

// Added by Scott J. Brandenberg
TclCharFn TclCommand_doPySimple1Gen;
TclCharFn TclCommand_doTzSimple1Gen;

// End added by SJB
// Added by Prishati Raychowdhury (UCSD)
TclCharFn TclSafeBuilder_doShallowFoundationGen;
// End PRC
//Leo
TclCharFn TclSafeBuilder_addRemoHFiber;
TclCharFn TclCommand_addFrictionModel;
TclCharFn TclCommand_addStiffnessDegradation;
TclCharFn TclCommand_addUnloadingRule;
TclCharFn TclCommand_addStrengthDegradation;
/// added by ZHY
TclCharFn TclCommand_UpdateMaterialStage;
TclCharFn TclCommand_UpdateMaterials;
/// added by ZHY
TclCharFn TclCommand_UpdateParameter;
////////////////gnp adding rayleigh /////////////////////
TclCharFn TclCommand_addElementRayleigh;
/////////////////////////////////////////////////////////

// Added by Alborz Ghofrani - U.Washington
TclCharFn TclCommand_GenerateInterfacePoints;

