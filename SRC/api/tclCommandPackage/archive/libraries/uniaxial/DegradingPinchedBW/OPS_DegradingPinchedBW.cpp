

#include <g3_api.h>


#include <SRC/material/uniaxial/DegradingPinchedBW.h>
void *OPS_DegradingPinchedBW(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData1[1];
  double dData[18];
  int iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DegradingPinchedBW tag"
           << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new DegradingPinchedBW(
      iData1[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
      dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
      dData[13], dData[14], dData[15], dData[16], dData[17], iData2[0]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "DegradingPinchedBW\n";
    return 0;
  }

  return theMaterial;
}
