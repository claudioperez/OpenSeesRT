

#include <g3_api.h>


#include <SRC/material/uniaxial/RambergOsgoodSteel.h>
void *OPS_RambergOsgoodSteel(void)
{
  if (numRambergOsgoodSteel == 0) {
    opserr
        << "RambergOsgoodSteel unaxial material - Written by R.Rahimi & "
           "R.Sepasdar & Dr. Mo. R. Banan Shiraz University Copyright 2012; \n";
    numRambergOsgoodSteel++;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[4];
  double mData[2];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial RambergOsgoodSteel tag"
           << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;
  }

  theMaterial =
      new RambergOsgoodSteel(iData[0], dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "RambergOsgoodSteel\n";
    return 0;
  }
  return theMaterial;
}
