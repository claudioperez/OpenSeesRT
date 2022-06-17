
#ifdef _RELIABILITY
// AddingSensitivity:BEGIN /////////////////////////////////////////////////
#  include <ReliabilityDomain.h>
#  include <SensitivityAlgorithm.h>
// AddingSensitivity:END /////////////////////////////////////////////////
#  include <TclReliabilityBuilder.h>
   int reliability(ClientData, Tcl_Interp *, int, TCL_Char **);
   int wipeReliability(ClientData, Tcl_Interp *, int, TCL_Char **);
   int optimization(ClientData, Tcl_Interp *, int, TCL_Char **); // Quan  (2)
#endif

#ifdef _RELIABILITY
// AddingSensitivity:BEGIN /////////////////////////////////////////////
static TclReliabilityBuilder *theReliabilityBuilder = 0;

Integrator *theSensitivityAlgorithm = 0;
Integrator *theSensitivityIntegrator = 0;
#include <TclOptimizationBuilder.h>
static TclOptimizationBuilder *theOptimizationBuilder =
    0; // Quan March 2010 (3)

// AddingSensitivity:END ///////////////////////////////////////////////
#endif



#ifdef _RELIABILITY
  Tcl_CreateCommand(interp, "wipeReliability", wipeReliability,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "reliability", reliability, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
  theReliabilityBuilder = 0;
  // AddingSensitivity:BEGIN //////////////////////////////////
  Tcl_CreateCommand(interp, "computeGradients", &computeGradients,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensitivityAlgorithm", &sensitivityAlgorithm,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensitivityIntegrator", &sensitivityIntegrator,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodeDisp", &sensNodeDisp, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensLambda", &sensLambda, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL); // Abbas
  Tcl_CreateCommand(interp, "sensNodeVel", &sensNodeVel, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodeAccel", &sensNodeAccel, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensSectionForce", &sensSectionForce,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodePressure", &sensNodePressure,
                    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  theSensitivityAlgorithm = 0;
  theSensitivityIntegrator = 0;
  // FMK RELIABILITY theReliabilityStaticAnalysis =0;
  // FMK RELIABILITY theReliabilityTransientAnalysis =0;
  // AddingSensitivity:END //////////////////////////////////

  theOptimizationBuilder = 0;

  // --- Quan March 2010  (4)
  Tcl_CreateCommand(interp, "optimization", &optimization, (ClientData)NULL,
                    (Tcl_CmdDeleteProc *)NULL);
#endif

#ifdef _RELIABILITY

// -- optimization Quan March 2010  (5)
int
optimization(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv)
{

  if (theOptimizationBuilder == 0) {

    theOptimizationBuilder = new TclOptimizationBuilder(theDomain, interp);

    return TCL_OK;
  } else
    return TCL_ERROR;
}

int
reliability(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  if (theReliabilityBuilder == 0) {

    theReliabilityBuilder = new TclReliabilityBuilder(theDomain, interp);
    return TCL_OK;
  } else
    return TCL_ERROR;
}

int
wipeReliability(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char **argv)
{
  if (theReliabilityBuilder != 0) {
    delete theReliabilityBuilder;
    theReliabilityBuilder = 0;
  }
  return TCL_OK;
}
// AddingSensitivity:END /////////////////////////////////////////////////

#endif


