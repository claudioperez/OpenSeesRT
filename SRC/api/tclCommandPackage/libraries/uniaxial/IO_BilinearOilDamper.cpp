

#include <g3_api.h>


#include <SRC/material/uniaxial/BilinearOilDamper.h>
void *OPS_BilinearOilDamper(void)
{
  if (numBilinearOilDamperMaterials == 0) {
    numBilinearOilDamperMaterials++;
    opserr << "BilinearOilDamper Model by Sarven Akcelyan and Dimitrios G. "
              "Lignos, PhD, McGill University\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[9];
  int numData = 1;
  // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  BilinearOilDamper tag"
           << endln;
    return 0;
  }
  // Check if the input variables
  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 2 && numData != 4 && numData != 5 && numData != 9) {
    opserr << "Invalid #args, want: uniaxialMaterial BilinearOilDamper "
           << iData[0]
           << " K? C? <Fr? p?> <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args want: uniaxialMaterial BilinearOilDamper "
           << iData[0]
           << " K? C? <Fr? p?> <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;

    return 0;
  }

  if (numData == 2) {
    // Default variables
    dData[2] = 1.0;
    dData[3] = 1.0;
    dData[4] = 0.0;
    dData[5] = 1;
    dData[6] = 0.000001;
    dData[7] = 0.0000000001;
    dData[8] = 15;
  }
  if (numData == 4) {
    // Default variables
    dData[4] = 0.0;
    dData[5] = 1;
    dData[6] = 0.000001;
    dData[7] = 0.0000000001;
    dData[8] = 15;
  }
  if (numData == 5) {
    // Default variables
    dData[5] = 1;
    dData[6] = 0.000001;
    dData[7] = 0.0000000001;
    dData[8] = 15;
  }

  // Parsing was successful, allocate the material with zero index
  theMaterial =
      new BilinearOilDamper(iData[0], dData[0], dData[1], dData[2], dData[3],
                            dData[4], dData[5], dData[6], dData[7], dData[8]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "BilinearOilDamper Material\n";
    return 0;
  }

  return theMaterial;
}
