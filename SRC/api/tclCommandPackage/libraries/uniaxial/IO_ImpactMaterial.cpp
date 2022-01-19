

#include <g3_api.h>


#include <SRC/material/uniaxial/ImpactMaterial.h>
void *OPS_ImpactMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();

  if (argc < 5) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial "
              "ImpactMaterial ?tag $K1 $K2 $Delta_y $gap"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[4];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ImpactMaterial tag" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: for ImpactMaterial tag: "
           << iData[0] << "\n";
    return 0;
  }

  theMaterial =
      new ImpactMaterial(iData[0], dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type ImpactMaterial\n";
    return 0;
  }

  return theMaterial;
}
