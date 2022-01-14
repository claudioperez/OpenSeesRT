

#include <g3_api.h>


#include <SRC/material/uniaxial/BWBN.h>
void *OPS_BWBN(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData1[1];
  double dData[13];
  int iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial BWBN tag" << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new BWBN(iData1[0], dData[0], dData[1], dData[2], dData[3],
                         dData[4], dData[5], dData[6], dData[7], dData[8],
                         dData[9], dData[10], dData[11], dData[12], iData2[0]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BWBN\n";
    return 0;
  }

  return theMaterial;
}
