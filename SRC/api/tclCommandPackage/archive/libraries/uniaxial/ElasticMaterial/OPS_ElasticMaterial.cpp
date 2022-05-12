

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticMaterial.h>
void *OPS_ElasticMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? <eta?> "
              "<Eneg?> ... "
           << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData >= 3) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
    }
  } else if (numData >= 2) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
    }
    dData[2] = dData[0];
  } else {
    numData = 1;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial Elastic " << iData[0]
             << endln;
      return 0;
    }
    dData[1] = 0.0;
    dData[2] = dData[0];
  }

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMaterial(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ElasticMaterial\n";
    return 0;
  }

  return theMaterial;
}
