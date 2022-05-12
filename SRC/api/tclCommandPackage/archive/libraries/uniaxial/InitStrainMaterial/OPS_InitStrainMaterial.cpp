

#include <g3_api.h>


#include <SRC/material/uniaxial/InitStrainMaterial.h>
OPS_Export void *OPS_InitStrainMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;

  int iData[2];
  double dData[1];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr
        << "WARNING invalid uniaxialMaterial InitStrainMaterial $tag $otherTag"
        << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "Could not find material with tag: " << iData[1]
           << "uniaxialMaterial InitStrain $tag $otherTag $eps0" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr
        << "Invalid Args want: uniaxialMaterial InitStrain $tag $otherTag $eps0"
        << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new InitStrainMaterial(iData[0], *theOtherMaterial, dData[0]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "InitStrainMaterial\n";
    return 0;
  }

  return theMaterial;
}
