

#include <g3_api.h>


#include <SRC/material/uniaxial/HookGap.h>
void *OPS_HookGap(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? gap? ... "
           << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial HookGapMaterial"
           << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData >= 3) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial HookGap " << iData[0]
             << endln;
      return 0;
    }
  } else {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial HookGap " << iData[0]
             << endln;
      return 0;
    }
    dData[1] = -dData[2];
    dData[2] = dData[1];
    ;
  }

  // Parsing was successful, allocate the material
  theMaterial = new HookGap(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type HookGap\n";
    return 0;
  }

  return theMaterial;
}
