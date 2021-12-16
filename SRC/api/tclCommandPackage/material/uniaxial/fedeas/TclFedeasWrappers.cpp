#include <tcl.h>
#include <stdio.h>
#include <elementAPI_G3.h>
#include <DegradingUniaxialWrapper.hh>

#define WRAPPER_CMD "FedeasDegradingWrapper"

int
TclSafeBuilder_addFedeasWrapper(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char **argv)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theWrappedMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain =  1.0e16;
  int tags[2];
  opserr << "hi";
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-min $minStrain> <-max $maxStrain>"
           << endln;
    return 0;
  }

  // Get wrapper tag
  if (Tcl_GetInt(interp, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return 0;
  }
  // Get base tag
  if (Tcl_GetInt(interp, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return 0;
  }

  // Get base material
  theWrappedMaterial = G3_getUniaxialMaterialInstance(interp, tags[1]);
  if (theWrappedMaterial == 0) {
    opserr << "WARNING invalid baseTag uniaxialMaterial " WRAPPER_CMD " tag: "
           << tags[1] << endln;
    return 0;
  }

  int argn = 2;
  while (argn < argc) {
    const char *param = argv[argn];
    double doubdat;

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-dmg") == 0)    ||
        (strcmp(param, "-DMG") == 0))   {
      
      opserr << "WARNING invalid baseTag uniaxialMaterial " WRAPPER_CMD ;
    } else {
      opserr << "WARNING invalid option:" << param
             << " uniaxialMaterial " WRAPPER_CMD " tag: " << endln;
      return 0;
    }
    argn++;
  }

// Parsing was successful, allocate the material
   theMaterial = new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial,
                                              minStrain, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type " WRAPPER_CMD << endln;
    return 0;
  }

  return G3_addUniaxialMaterial(interp, theMaterial);
  // return theMaterial;
}

