

#include <g3_api.h>


#include <SRC/material/uniaxial/Maxwell.h>
void *OPS_Maxwell(G3_Runtime *rt)
{
  if (numMaxwellMaterials == 0) {
    numMaxwellMaterials++;
    opserr << "Maxwell Model - D.Lignos, McGill University\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[2];
  double dData[4];
  iData[1] = 0;

  int numRData = OPS_GetNumRemainingInputArgs();
  if (numRData != 5 && numRData != 6) {
    opserr << "Invalid #args for command uniaxialMaterial Maxwell\n";
    return 0;
  }

  // Check tag
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Maxwell tag" << endln;
    return 0;
  }
  // Check if we have 4 input variables for K, C, Alpha, L
  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial Maxwell tag? K? C? Alpha? "
              "Length L?"
           << endln;
    return 0;
  }

  if (numRData == 6) {
    const char *cArray = OPS_GetString();
    // OPS_GetStringCopy(&cArray);
    if ((strcmp(cArray, "-returnD") == 0) || (strcmp(cArray, "-D") == 0))
      iData[1] = 1;
    delete[] cArray;
  }

  // Parsing was successful, allocate the material with zero index
  theMaterial =
      new Maxwell(iData[0], dData[0], dData[1], dData[2], dData[3], iData[1]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Maxwell "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
