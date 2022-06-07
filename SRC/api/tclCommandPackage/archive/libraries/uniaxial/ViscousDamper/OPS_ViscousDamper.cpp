

#include <g3_api.h>


#include <SRC/material/uniaxial/ViscousDamper.h>
void *OPS_ViscousDamper(void)

{
  if (numViscousDamperMaterials == 0) {
    numViscousDamperMaterials++;
    opserr << "ViscousDamper Model by Sarven Akcelyan and Dimitrios G. Lignos, "
              "PhD, McGill University\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[8];
  int numData = 1;
  // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ViscousDamper tag" << endln;
    return 0;
  }
  // Check if we have 3 or 4 or 8 input variables
  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 4 && numData != 8) {
    opserr << "Invalid #args, want: uniaxialMaterial ViscousDamper " << iData[0]
           << " K? C? Alpha? <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args want: uniaxialMaterial ViscousDamper " << iData[0]
           << " K? C? Alpha? <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;

    return 0;
  }

  if (numData == 3) {
    // Default variables
    dData[3] = 0.0;
    dData[4] = 1;
    dData[5] = 0.000001;
    dData[6] = 0.0000000001;
    dData[7] = 15;
  }

  if (numData == 4) {
    // Default variables
    dData[4] = 1;
    dData[5] = 0.000001;
    dData[6] = 0.0000000001;
    dData[7] = 15;
  }

  // Parsing was successful, allocate the material with zero index
  theMaterial =
      new ViscousDamper(iData[0], dData[0], dData[1], dData[2], dData[3],
                        dData[4], dData[5], dData[6], dData[7]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ViscousDamper "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
