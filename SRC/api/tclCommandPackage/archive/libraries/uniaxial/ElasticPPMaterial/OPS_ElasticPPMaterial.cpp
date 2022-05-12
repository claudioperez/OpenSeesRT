

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticPPMaterial.h>
void *OPS_ElasticPPMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3 || numArgs > 5) {
    opserr << "Invalid #args,  want: uniaxialMaterial ElasticPP $tag $E $epsP "
              "<$epsN $eps0>\n";
    return 0;
  }

  int iData[1];
  double dData[4];
  dData[3] = 0.0; // setting default eps0 to 0.

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial ElasticPP" << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial ElasticPP " << iData[0] << endln;
    return 0;
  }

  if (numData == 2)
    dData[2] = -dData[1];

  // Parsing was successful, allocate the material
  theMaterial =
      new ElasticPPMaterial(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticPP\n";
    return 0;
  }

  return theMaterial;
}
