

#include <g3_api.h>


#include <SRC/material/uniaxial/SeriesMaterial.h>
void *OPS_SeriesMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr
        << "Invalid #args,  want: uniaxialMaterial Series $tag $tag1 $tag2 ... "
        << endln;
    return 0;
  }

  int *iData = new int[numArgs];
  UniaxialMaterial **theMats = new UniaxialMaterial *[numArgs - 1];

  if (OPS_GetIntInput(&numArgs, iData) != 0) {
    opserr << "WARNING invalid data for uniaxialMaterial Series" << endln;
    return 0;
  }

  for (int i = 1; i < numArgs; i++) {
    UniaxialMaterial *theMat = OPS_GetUniaxialMaterial(iData[i]);
    if (theMat == 0) {
      opserr << "WARNING no existing material with tag " << iData[i]
             << " for uniaxialMaterial Series" << iData[0] << endln;
      delete[] iData;
      delete[] theMats;
      return 0;
    }
    theMats[i - 1] = theMat;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SeriesMaterial(iData[0], numArgs - 1, theMats);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Series\n";
    return 0;
  }

  delete[] iData;
  delete[] theMats;

  return theMaterial;
}
