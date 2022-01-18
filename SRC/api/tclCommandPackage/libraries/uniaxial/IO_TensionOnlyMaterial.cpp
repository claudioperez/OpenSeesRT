

#include <g3_api.h>


#include <SRC/material/uniaxial/TensionOnlyMaterial.h>
OPS_Export void *OPS_TensionOnlyMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial TensionOnly $tag $otherTag "
              "<-min $minStrain> <-max $maxStrain>"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial TensionOnly $tag $otherTag"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial TensionOnly tag: "
           << iData[0] << endln;
    return 0;
  }

  argc = OPS_GetNumRemainingInputArgs();
  while (argc > 1) {
    // char argvLoc[10];
    const char *argvLoc = OPS_GetString();
    /*    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid string option uniaxialMaterial TensionOnly tag:
    " << iData[0] << endln; return 0;
    }
    */
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) ||
        (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {
        opserr
            << "WARNING invalid min value  uniaxialMaterial TensionOnly tag: "
            << iData[0] << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) ||
               (strcmp(argvLoc, "-Max") == 0) ||
               (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {
        opserr
            << "WARNING invalid min value  uniaxialMaterial TensionOnly tag: "
            << iData[0] << endln;
        return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc
             << " uniaxialMaterial TensionOnly tag: " << iData[0] << endln;
      return 0;
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new TensionOnlyMaterial(iData[0], *theOtherMaterial);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "TensionOnlyMaterial\n";
    return 0;
  }

  return theMaterial;
}
