

#include <g3_api.h>


#include <SRC/material/uniaxial/SimpleFractureMaterial.h>
void *OPS_SimpleFractureMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double maxStrain = 1.0e16;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3) {
    opserr << "WARNING invalid uniaxialMaterial SimpleFracture $tag $otherTag "
              "$maxStrain>"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SimpleFracture $tag $otherTag "
              "$maxStrain"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag:  uniaxialMaterial SimpleFracture $tag "
              "$otherTag $max: "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &maxStrain) != 0) {
    opserr << "WARNING invalid maxStrain: uniaxialMaterial  SimpleFracture "
              "$tag $otherTag $maxStrain"
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new SimpleFractureMaterial(iData[0], *theOtherMaterial, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "SimpleFractureMaterial\n";
    return 0;
  }

  return theMaterial;
}
