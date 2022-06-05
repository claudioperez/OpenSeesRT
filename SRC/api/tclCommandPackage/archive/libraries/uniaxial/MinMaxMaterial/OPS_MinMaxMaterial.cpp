

#include <g3_api.h>


#include <SRC/material/uniaxial/MinMaxMaterial.h>
OPS_Export void *OPS_MinMaxMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial MinMaxMaterial $tag $otherTag "
              "<-min $minStrain> <-max $maxStrain>"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MinMaxMaterial $tag $otherTag"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial MinMax tag: "
           << iData[0] << endln;
    return 0;
  }

  argc = OPS_GetNumRemainingInputArgs();
  while (argc > 1) {
    // char argvLoc[10];
    const char *argvLoc = OPS_GetString();
    /*    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid string option uniaxialMaterial MinMax tag: " <<
    iData[0] << endln; return 0;
    }
    */
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) ||
        (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {
        opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: "
               << iData[0] << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) ||
               (strcmp(argvLoc, "-Max") == 0) ||
               (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {
        opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: "
               << iData[0] << endln;
        return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc
             << " uniaxialMaterial MinMax tag: " << iData[0] << endln;
      return 0;
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new MinMaxMaterial(iData[0], *theOtherMaterial, minStrain, maxStrain);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type MinMaxMaterial\n";
    return 0;
  }

  return theMaterial;
}
