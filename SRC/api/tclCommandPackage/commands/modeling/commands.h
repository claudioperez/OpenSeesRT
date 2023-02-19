/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// cmp
//
#include <tcl.h>

// modeling/model.cpp
extern Tcl_CmdProc  TclCommand_wipeModel;
extern Tcl_CmdProc  buildModel;

// modeling/nodes.cpp
extern Tcl_CmdProc  TclCommand_getNDM;
extern Tcl_CmdProc  TclCommand_getNDF;
extern Tcl_CmdProc  TclCommand_addNode;
extern Tcl_CmdProc  TclCommand_addNodalMass;
// 
extern Tcl_CmdProc  TclCommand_addPattern;
extern Tcl_CmdProc  TclCommand_addTimeSeries;
extern Tcl_CmdProc  TclCommand_addGeomTransf;

extern Tcl_CmdProc  TclCommand_addElement;
extern Tcl_CmdProc  TclCommand_doBlock2D;
extern Tcl_CmdProc  TclCommand_doBlock3D;
extern Tcl_CmdProc  TclCommand_addUniaxialMaterial;
extern Tcl_CmdProc  TclCommand_addHystereticBackbone;
extern Tcl_CmdProc  TclCommand_addSection;
extern Tcl_CmdProc  TclCommand_addPatch;
extern Tcl_CmdProc  TclCommand_addReinfLayer;
// extern Tcl_CmdProc  TclCommand_addRemoFiber;
extern Tcl_CmdProc  TclCommand_addFiber;
// Constraints
extern Tcl_CmdProc  TclCommand_addHomogeneousBC;
extern Tcl_CmdProc  TclCommand_addEqualDOF_MP;
// Loads
// extern Tcl_CmdProc  TclCommand_addNodalLoad;
Tcl_CmdProc TclCommand_addElementalLoad;

// Damping
Tcl_CmdProc modalDamping;
Tcl_CmdProc modalDampingQ;

// Constraints
Tcl_CmdProc TclCommand_addHomogeneousBC_X;
Tcl_CmdProc TclCommand_addHomogeneousBC_Y; 
Tcl_CmdProc TclCommand_addHomogeneousBC_Z;
Tcl_CmdProc TclCommand_addEqualDOF_MP_Mixed;
Tcl_CmdProc TclCommand_addMP;
Tcl_CmdProc TclCommand_addSP;
Tcl_CmdProc TclCommand_RigidLink;
Tcl_CmdProc TclCommand_addImposedMotionSP;
Tcl_CmdProc TclCommand_addGroundMotion;
Tcl_CmdProc TclCommand_RigidDiaphragm;

// Other
Tcl_CmdProc TclCommand_addParameter;
Tcl_CmdProc TclCommand_mesh;
Tcl_CmdProc TclCommand_remesh;
Tcl_CmdProc TclCommand_backgroundMesh; 
Tcl_CmdProc TclCommand_addBeamIntegration;

Tcl_CmdProc TclCommand_addLimitCurve;
Tcl_CmdProc TclCommand_addNDMaterial;
Tcl_CmdProc TclCommand_addSeries;


struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
  bool was_added = false;

}  const tcl_char_cmds[] =  {
  {"wipe",             TclCommand_wipeModel},
  {"build",            buildModel},

  {"getNDM",           TclCommand_getNDM},
  {"getNDF",           TclCommand_getNDF},
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

//   {"load",             TclCommand_addNodalLoad},
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

  {"eleLoad",          TclCommand_addElementalLoad},

//{"beamIntegration",  TclCommand_addBeamIntegration},

/*
  {"mp",                   TclCommand_addMP},

  {"block2D",              TclCommand_doBlock2D},
  {"block3D",              TclCommand_doBlock3D},

  {"equalDOF_Mixed",       TclCommand_addEqualDOF_MP_Mixed},
  {"rigidDiaphragm",       &TclCommand_RigidDiaphragm},
  {"PySimple1Gen",         TclCommand_doPySimple1Gen},
  {"TzSimple1Gen",         TclCommand_doTzSimple1Gen},
  {"ShallowFoundationGen", BasicModelBuilder_doShallowFoundationGen},
  {"Hfiber",               BasicModelBuilder_addRemoHFiber},
*/

// OTHER OBJECT TYPES
  {"hystereticBackbone",   TclCommand_addHystereticBackbone},
  {          "backbone",   TclCommand_addHystereticBackbone},

#if 0
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
#endif

#if 0
  {"parameter",            TclCommand_addParameter},
  {"addToParameter",       TclCommand_addParameter},
  {"updateParameter",      TclCommand_addParameter},
#endif

#if 0
// command for elast2plast in Multi-yield plasticity, by ZHY
  {"updateMaterialStage", TclCommand_UpdateMaterialStage},
  {"updateMaterials",     TclCommand_UpdateMaterials},
#endif

#if 0
// command for updating properties of soil materials, by ZHY
   {"updateParameter", TclCommand_UpdateParameter},
#endif

};



Tcl_CmdProc TclCommand_Package;

// Added by Scott J. Brandenberg
Tcl_CmdProc TclCommand_doPySimple1Gen;
Tcl_CmdProc TclCommand_doTzSimple1Gen;

// End added by SJB
// Added by Prishati Raychowdhury (UCSD)
Tcl_CmdProc BasicModelBuilder_doShallowFoundationGen;
// End PRC
//Leo
Tcl_CmdProc BasicModelBuilder_addRemoHFiber;
Tcl_CmdProc TclCommand_addFrictionModel;
Tcl_CmdProc TclCommand_addStiffnessDegradation;
Tcl_CmdProc TclCommand_addUnloadingRule;
Tcl_CmdProc TclCommand_addStrengthDegradation;
/// added by ZHY
Tcl_CmdProc TclCommand_UpdateMaterialStage;
Tcl_CmdProc TclCommand_UpdateMaterials;
/// added by ZHY
Tcl_CmdProc TclCommand_UpdateParameter;
////////////////gnp adding rayleigh /////////////////////
Tcl_CmdProc TclCommand_addElementRayleigh;
/////////////////////////////////////////////////////////

// Added by Alborz Ghofrani - U.Washington
Tcl_CmdProc TclCommand_GenerateInterfacePoints;

