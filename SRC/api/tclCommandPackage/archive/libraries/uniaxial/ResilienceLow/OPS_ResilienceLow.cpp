

#include <g3_api.h>


#include <SRC/material/uniaxial/ResilienceLow.h>
void *OPS_ResilienceLow(void)
{

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int iData[1];
  double dData[5];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterialtag" << endln;
    return 0;
  }
  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 5) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceLow " << iData[0]
           << "  PY DPmax Pmax Ke Kd" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceLow " << iData[0]
           << "  PY DPmax Pmax Ke Kd" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new ResilienceLow(iData[0], dData[0], dData[1], dData[2],
                                  dData[3], dData[4]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type ResilienceLow\n";
    return 0;
  }

  return theMaterial;
}
