/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <runtimeAPI.h>
class G3_Runtime;

OPS_Routine OPS_NewtonRaphsonAlgorithm;
OPS_Routine OPS_ExpressNewton;
OPS_Routine OPS_ModifiedNewton;
OPS_Routine OPS_NewtonHallM;

OPS_Routine OPS_Newmark;
OPS_Routine OPS_StagedNewmark;
OPS_Routine OPS_GimmeMCK;
OPS_Routine OPS_AlphaOS;
OPS_Routine OPS_AlphaOS_TP;
OPS_Routine OPS_AlphaOSGeneralized;
OPS_Routine OPS_AlphaOSGeneralized_TP;
OPS_Routine OPS_ExplicitDifference;
OPS_Routine OPS_CentralDifference;
OPS_Routine OPS_CentralDifferenceAlternative;
OPS_Routine OPS_CentralDifferenceNoDamping;
OPS_Routine OPS_Collocation;
OPS_Routine OPS_CollocationHSFixedNumIter;
OPS_Routine OPS_CollocationHSIncrLimit;
OPS_Routine OPS_CollocationHSIncrReduct;
OPS_Routine OPS_GeneralizedAlpha;
OPS_Routine OPS_HHT;
OPS_Routine OPS_HHT_TP;
OPS_Routine OPS_HHTExplicit;
OPS_Routine OPS_HHTExplicit_TP;
OPS_Routine OPS_HHTGeneralized;
OPS_Routine OPS_HHTGeneralized_TP;
OPS_Routine OPS_HHTGeneralizedExplicit;
OPS_Routine OPS_HHTGeneralizedExplicit_TP;
OPS_Routine OPS_HHTHSFixedNumIter;
OPS_Routine OPS_HHTHSFixedNumIter_TP;
OPS_Routine OPS_HHTHSIncrLimit;
OPS_Routine OPS_HHTHSIncrLimit_TP;
OPS_Routine OPS_HHTHSIncrReduct;
OPS_Routine OPS_HHTHSIncrReduct_TP;
OPS_Routine OPS_KRAlphaExplicit;
OPS_Routine OPS_KRAlphaExplicit_TP;
OPS_Routine OPS_NewmarkExplicit;
OPS_Routine OPS_NewmarkHSFixedNumIter;
OPS_Routine OPS_NewmarkHSIncrLimit;
OPS_Routine OPS_NewmarkHSIncrReduct;
OPS_Routine OPS_WilsonTheta;

// wipe is not added to the table on purpose
static Tcl_CmdProc wipeAnalysis;
static Tcl_CmdProc specifyAnalysis;
static Tcl_CmdProc eigenAnalysis;
static Tcl_CmdProc modalProperties;
static Tcl_CmdProc responseSpectrum;
static Tcl_CmdProc printA;
static Tcl_CmdProc printB;
static Tcl_CmdProc initializeAnalysis;
static Tcl_CmdProc resetModel;
static Tcl_CmdProc analyzeModel;
static Tcl_CmdProc specifyConstraintHandler;
// Damping
static Tcl_CmdProc modalDamping;

extern Tcl_CmdProc specifyIntegrator;

extern Tcl_CmdProc specifySOE;
extern Tcl_CmdProc specifySysOfEqnTable;

// commands/analysis/algorithm.cpp
extern Tcl_CmdProc TclCommand_specifyAlgorithm;
extern Tcl_CmdProc TclCommand_numIter;
extern Tcl_CmdProc TclCommand_accelCPU;
extern Tcl_CmdProc TclCommand_totalCPU;
extern Tcl_CmdProc TclCommand_solveCPU;
extern Tcl_CmdProc TclCommand_numFact;

// from commands/analysis/ctest.cpp
extern Tcl_CmdProc specifyCTest;
extern Tcl_CmdProc getCTestNorms;
extern Tcl_CmdProc getCTestIter;
extern Tcl_CmdProc TclCommand_algorithmRecorder;

struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
}  const tcl_analysis_cmds[] =  {
    {"system",            &specifySysOfEqnTable},

    {"test",              &specifyCTest},
    {"testIter",          &getCTestIter},
    {"testNorms",         &getCTestNorms},
    {"integrator",        &specifyIntegrator},
    {"constraints",       &specifyConstraintHandler},

    {"eigen",             &eigenAnalysis},
    {"analysis",          &specifyAnalysis},

    {"analyze",           &analyzeModel},
    {"initialize",        &initializeAnalysis},
    {"modalProperties",   &modalProperties},
    {"modalDamping",      &modalDamping},
    {"modalDampingQ",     &modalDamping},
    {"responseSpectrum",  &responseSpectrum},
    {"printA",            &printA},
    {"printB",            &printB},
    {"reset",             &resetModel},

  // From algorithm.cpp
    {"algorithm", &TclCommand_specifyAlgorithm},
    {"numIter",   &TclCommand_numIter},
    {"numFact",   &TclCommand_numFact},
    {"accelCPU",  &TclCommand_accelCPU},
    {"totalCPU",  &TclCommand_totalCPU},
    {"solveCPU",  &TclCommand_solveCPU},
  // recorder.cpp
    {"algorithmRecorder",   &TclCommand_algorithmRecorder},
};

