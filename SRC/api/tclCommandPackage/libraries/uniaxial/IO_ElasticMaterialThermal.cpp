

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticMaterialThermal.h>
void *OPS_ElasticMaterialThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int softindex = 0;
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? "
              "alpha?<eta?> ... "
           << endln;
    return 0;
  }

  int iData[1];
  double dData1[2];
  double dData2[2];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData == 1)
    dData1[1] = 0.0;
  else
    numData = 2;

  if (OPS_GetDoubleInput(&numData, dData1) != 0) {
    opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData > 0) {
    const char *typeChar = OPS_GetString();
    if ((strcmp(typeChar, "-SteelSoft") == 0) ||
        (strcmp(typeChar, "-SSoft") == 0) ||
        (strcmp(typeChar, "-sSoft") == 0)) {
      softindex = 1;
    } else if ((strcmp(typeChar, "-ConcreteSoft") == 0) ||
               (strcmp(typeChar, "-CSoft") == 0) ||
               (strcmp(typeChar, "-cSoft") == 0)) {
      softindex = 2;
    }
  }

  dData2[0] = 0.0;
  dData2[1] = 0.0;

  numData = numData - 1;
  if (numData > 2)
    numData = 2;
  if (numData > 0) {
    if (OPS_GetDoubleInput(&numData, dData2) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
    }
  }
  // opserr<< "received alpha"<<dData[1]<<endln;

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMaterialThermal(iData[0], dData1[0], dData1[1],
                                           dData2[0], dData2[1], softindex);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ElasticMaterialThermal\n";
    return 0;
  }

  return theMaterial;
}
