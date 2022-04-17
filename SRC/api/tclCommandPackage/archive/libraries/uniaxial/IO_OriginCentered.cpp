

#include <g3_api.h>


#include <SRC/material/uniaxial/OriginCentered.h>
void *OPS_OriginCentered(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial OriginCentered tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 6) {
    opserr << "Invalid #args, want: uniaxialMaterial OriginCentered "
           << iData[0] << " f1? e1? f2? e2? f3? e3?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid arggs: uniaxialMaterial OriginCentered " << iData[0]
           << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new OriginCentered(iData[0], dData[0], dData[1], dData[2],
                                   dData[3], dData[4], dData[5]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "OriginCentered Material\n";
    return 0;
  }

  return theMaterial;
}
