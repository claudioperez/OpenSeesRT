

#include <g3_api.h>


#include <SRC/material/uniaxial/BraceMaterial.h>
void *OPS_BraceMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[17];
  int numData = 1;
  // Check Tag and number of Fingers
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Cast Fuse tag" << endln;
    return 0;
  }

  int numRemaining = OPS_GetNumRemainingInputArgs();
  if (numRemaining != 17 || numRemaining != 13) {
    if (OPS_GetDoubleInput(&numRemaining, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial CastFuse tag? NLegs? bo? "
                "h? Fy? E? L? b? R0? cR1? cR2? a1? a2? a3? a4?";
      return 0;
    }

    // Parsing was successful, allocate the material
    if (numRemaining == 13) {
      theMaterial = new BraceMaterial(iData[0], dData[0], dData[1], dData[2],
                                      dData[3], dData[4], dData[5], dData[6],
                                      dData[7], dData[8], dData[9], dData[10],
                                      dData[11], dData[12], dData[12]);
    } else {
      theMaterial = new BraceMaterial(
          iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
          dData[6], dData[7], dData[8], dData[9], dData[10], dData[11],
          dData[12], dData[13], dData[14], dData[15], dData[16]);
    }
  }

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type Cast Material\n";
    return 0;
  }

  return theMaterial;
}
