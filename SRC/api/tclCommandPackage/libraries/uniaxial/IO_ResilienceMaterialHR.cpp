

#include <g3_api.h>


#include <SRC/material/uniaxial/ResilienceMaterialHR.h>
void *OPS_ResilienceMaterialHR(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int iData[1];
  double dData[7];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ResilienceMaterialHR tag"
           << endln;
    return 0;
  }
  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceMaterialHR "
           << iData[0] << " DY PY DPmax Pmax Ke Kd coefficient" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceMaterialHR "
           << iData[0] << " DY PY DPmax Pmax Ke Kd coefficient" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new ResilienceMaterialHR(iData[0], dData[0], dData[1], dData[2], dData[3],
                               dData[4], dData[5], dData[6]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ResilienceMaterialHR\n";
    return 0;
  }

  return theMaterial;
}
