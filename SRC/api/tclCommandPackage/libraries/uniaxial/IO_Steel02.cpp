

#include <g3_api.h>
#include <SRC/material/uniaxial/Steel02.h>

void *OPS_Steel02(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 6 && numData != 10 && numData != 11) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02 " << iData[0]
           << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid double: uniaxialMaterial Steel02 " << iData[0]
             << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2]);

  } else if (numData == 6) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid int: uniaxialMaterial Steel02 " << iData[0]
             << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2], dData[3],
                              dData[4], dData[5]);

  } else if (numData == 10) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0]
             << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial =
        new Steel02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                    dData[5], dData[6], dData[7], dData[8], dData[9]);

  } else if (numData == 11) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0]
             << " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02(iData[0], dData[0], dData[1], dData[2], dData[3],
                              dData[4], dData[5], dData[6], dData[7], dData[8],
                              dData[9], dData[10]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel02 "
              "Material\n";
    return 0;
  }

  return theMaterial;
}

