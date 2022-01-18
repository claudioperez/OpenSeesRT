

#include <g3_api.h>


#include <SRC/material/uniaxial/SteelMPF.h>
OPS_Export void *OPS_SteelMPF(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  // Parse the script for material parameters
  if (numArgs != 9 && numArgs != 13) {
    opserr
        << "Incorrect # args, Want: uniaxialMaterial SteelMPF tag? sigyieldp? "
           "sigyieldn? E0? bp? bn? R0? cR1? cR2? <a1? a2? a3? a4?>";
    return 0;
  }

  int iData[1];
  double dData[12];
  dData[8] = 0.0;  // set a3 to constructor default
  dData[9] = 1.0;  // set a4 to constructor default
  dData[10] = 0.0; // set a5 to constructor default
  dData[11] = 1.0; // set a6 to constructor default

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelMPF tag" << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxialMaterial SteelMPF " << dData[0]
           << endln;
    return 0;
  }

  theMaterial = new SteelMPF(iData[0], dData[0], dData[1], dData[2], dData[3],
                             dData[4], dData[5], dData[6], dData[7], dData[8],
                             dData[9], dData[10], dData[11]);

  return theMaterial;
}
