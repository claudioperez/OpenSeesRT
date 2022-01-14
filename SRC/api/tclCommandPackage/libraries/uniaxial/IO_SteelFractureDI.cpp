

#include <g3_api.h>


#include <SRC/material/uniaxial/SteelFractureDI.h>
void *OPS_SteelFractureDI(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[14];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelFractureDI tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 14) {
    opserr
        << "Invalid #args, want: uniaxialMaterial SteelFractureDI " << iData[0]
        << " fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m? sigmin? FI_lim?"
        << endln;
    return 0;
  } else {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial SteelFractureDI " << iData[0]
             << " fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m? sigmin? "
                "FI_lim?"
             << endln;
      return 0;
    }
    theMaterial = new SteelFractureDI(iData[0], dData[0], dData[1], dData[2],
                                      dData[3], dData[4], dData[5], dData[6],
                                      dData[7], dData[8], dData[9], dData[10],
                                      dData[11], dData[12], dData[13]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "SteelFractureDI Material\n";
    return 0;
  }

  return theMaterial;
}
