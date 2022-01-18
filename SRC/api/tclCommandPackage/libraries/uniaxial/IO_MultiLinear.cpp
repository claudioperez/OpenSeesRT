

#include <g3_api.h>


#include <SRC/material/uniaxial/MultiLinear.h>
void *OPS_MultiLinear(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Invalid #args,  want: uniaxialMaterial MultiLinear tag? e1 s1 "
              "e2 s2 ... "
           << endln;
    return 0;
  }

  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag or soilType uniaxialMaterial "
              "MultiLinearMaterial"
           << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  int numSlope = numData / 2;
  double *dData = new double[numData];
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid pyData data for material uniaxial MultiLinear "
           << iData[0] << endln;
    return 0;
  }

  Vector e(numSlope);
  Vector s(numSlope);
  for (int i = 0; i < numSlope; i++) {
    e(i) = dData[2 * i];
    s(i) = dData[2 * i + 1];
  }

  // Parsing was successful, allocate the material
  theMaterial = new MultiLinear(iData[0], s, e);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MultiLinear\n";
    return 0;
  }

  return theMaterial;
}
