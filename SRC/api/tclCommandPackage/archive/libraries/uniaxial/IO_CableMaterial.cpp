

#include <g3_api.h>


#include <SRC/material/uniaxial/CableMaterial.h>
void *OPS_CableMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 5) {
    opserr << "Invalid # args, want: uniaxialMaterial Cable tag? $presetress "
              "$E $effUnitWeight $Lelement \n";
    return 0;
  }

  int iData[1];
  double dData[4];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Cable" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Cable " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new CableMaterial(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Cable\n";
    return 0;
  }

  return theMaterial;
}
