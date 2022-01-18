

#include <g3_api.h>


#include <SRC/material/uniaxial/Bond_SP01.h>
void *OPS_Bond_SP01(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numInput = OPS_GetNumRemainingInputArgs();
  if (numInput != 7 && numInput != 11) {
    opserr << "Invalid #args,  uniaxialMaterial Bond_SP01 tag? fy? sy? fu? su? "
              "b? R?";
    opserr << " <Cd? db? fc? la?>" << endln;
    return 0;
  }

  int iData[1];
  double dData[10];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = numInput - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  if (numInput == 7)
    theMaterial = new Bond_SP01(iData[0], dData[0], dData[1], dData[2],
                                dData[3], dData[4], dData[5]);
  else
    theMaterial = new Bond_SP01(iData[0], dData[0], dData[1], dData[2],
                                dData[3], dData[4], dData[5], dData[6],
                                dData[7], dData[8], dData[9]);

  return theMaterial;
}
