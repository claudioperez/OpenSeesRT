

#include <g3_api.h>


#include <SRC/material/uniaxial/MultiplierMaterial.h>
OPS_Export void *OPS_MultiplierMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING insufficient args, uniaxialMaterial Multiplier $tag "
              "$otherTag $multiplier"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Multiplier $tag $otherTag "
              "$multiplier"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial Multiplier tag: "
           << iData[0] << endln;
    return 0;
  }

  double multiplier = 1.0;
  numData = 1;
  if (OPS_GetDouble(&numData, &multiplier) < 0) {
    opserr << "WARNING invalid input uniaxialMaterial Multiplier tag: "
           << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new MultiplierMaterial(iData[0], *theOtherMaterial, multiplier);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "MultiplierMaterial\n";
    return 0;
  }

  return theMaterial;
}
