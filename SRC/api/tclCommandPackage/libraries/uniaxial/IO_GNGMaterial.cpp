

#include <g3_api.h>


#include <SRC/material/uniaxial/GNGMaterial.h>
void *OPS_GNGMaterial(G3_Runtime *rt)
{

  if (numGNGMaterials == 0) {
    numGNGMaterials++;
    opserr << "Grip 'n' Grab device installed in this structure!\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 4) {
    opserr << "Invalid #args,  want: uniaxialMaterial GNG tag E sigY P <eta>\n";
    return 0;
  }

  int tag;
  double dData[4];
  dData[3] = 0.0; // setting default eta to 0.

  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial GNG" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData > 4)
    numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial GNG \n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new GNGMaterial(tag, dData[0], dData[1], dData[2], dData[3]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type GNG\n";
    return 0;
  }

  return theMaterial;
}
