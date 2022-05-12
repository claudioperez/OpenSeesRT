#include <g3_api.h>
#include <SRC/material/uniaxial/Cast.h>

void *OPS_Cast(void)
{
  if (numCastMaterials == 0) {
    numCastMaterials++;
    opserr << "Cast Fuse uniaxial material - Written by Dimitrios G. Lignos, "
              "Ph.D.\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[14];
  int numData = 1;
  // Check Tag and number of Fingers
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Cast Fuse tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData < 14) {
    opserr << "WARNING insufficient number of args want  uniaxialMaterial "
              "CastFuse tag? NLegs? bo? h? Fy? E? L? b? R0? cR1? cR2? a1? a2? "
              "a3? a4\n";
    return 0;
  }
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial CastFuse tag? NLegs? bo? h? "
              "Fy? E? L? b? R0? cR1? cR2? a1? a2? a3? a4?";
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Cast(iData[0], dData[0], dData[1], dData[2], dData[3],
                         dData[4], dData[5], dData[6], dData[7], dData[8],
                         dData[9], dData[10], dData[11], dData[12], dData[13]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type Cast Material\n";
    return 0;
  }

  return theMaterial;
}
