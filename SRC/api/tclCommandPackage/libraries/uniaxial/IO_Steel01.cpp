

#include <g3_api.h>


#include <SRC/material/uniaxial/Steel01.h>
void *OPS_Steel01(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0]
           << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0]
           << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    dData[3] = STEEL_01_DEFAULT_A1;
    dData[4] = STEEL_01_DEFAULT_A2;
    dData[5] = STEEL_01_DEFAULT_A3;
    dData[6] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Steel01(iData[0], dData[0], dData[1], dData[2], dData[3],
                            dData[4], dData[5], dData[6]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel01 "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
