

#include <g3_api.h>


#include <SRC/material/uniaxial/Steel02Fatigue.h>
void *OPS_Steel02Fatigue(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[18];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02Fatigue tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 9 && numData != 12 && numData != 16 && numData != 17) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02Fatigue "
           << iData[0]
           << " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? cR1? "
              "cR2? <a1? a2? a3? a4?>>"
           << endln;
    return 0;
  }

  if (numData == 9) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid double: uniaxialMaterial Steel02Fatigue " << iData[0]
             << " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? "
                "cR1? cR2? <a1? a2? a3? a4?>>"
             << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial =
        new Steel02Fatigue(iData[0], dData[0], dData[1], dData[2], dData[3],
                           dData[4], dData[5], dData[6], dData[7], dData[8]);

  } else if (numData == 12) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid int: uniaxialMaterial Steel02Fatigue " << iData[0]
             << " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? "
                "cR1? cR2? <a1? a2? a3? a4?>>"
             << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);

  } else if (numData == 16) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Fatigue " << iData[0]
             << " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? "
                "cR1? cR2? <a1? a2? a3? a4?>>"
             << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15]);

  } else if (numData == 17) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02Fatigue " << iData[0]
             << " fy? E? b? Cd? Cf? alpha? beta? minStrain? maxStrain? <R0? "
                "cR1? cR2? <a1? a2? a3? a4?>>"
             << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Steel02Fatigue(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "Steel02Fatigue Material\n";
    return 0;
  }

  return theMaterial;
}
