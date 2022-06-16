#include <g3_api.h>
#include <tcl.h>
#include <Domain.h>
#include <UniaxialMaterial.h>
#include <WrapperUniaxialMaterial.h>

// extern TCL_Char **currentArgv;
// extern int currentArg;
// extern int maxArg;

UniaxialMaterial *
Tcl_addWrapperUniaxialMaterial(matObj *theMat, ClientData clientData,
                               Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // cmp //  theInterp = interp;

  // currentArgv = argv;
  // currentArg = 2;
  // maxArg = argc;

  G3_Runtime* rt = G3_getRuntime(interp);
  Domain* theDomain = G3_getDomain(rt);

  // get the current load factor
  static modelState theModelState;
  if (theDomain != 0) {
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;
    theModelState.time = time;
    theModelState.dt = dt;
  }

  // invoke the mat function with isw = 0
  int isw = ISW_INIT;
  int result = 0;
  theMat->matFunctPtr(theMat, &theModelState, 0, 0, 0, &isw, &result);
  int matType = theMat->matType; // GR added to support material

  if (result != 0 || matType != OPS_UNIAXIAL_MATERIAL_TYPE) {
    opserr << "Tcl_addWrapperUniaxialMaterial - failed in element function "
           << result << endln;
    return 0;
  }

  WrapperUniaxialMaterial *theMaterial =
      new WrapperUniaxialMaterial(argv[1], theMat);

  return theMaterial;
}
