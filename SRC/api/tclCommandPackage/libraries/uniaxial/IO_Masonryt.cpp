

#include <g3_api.h>


#include <SRC/material/uniaxial/Masonryt.h>
void *OPS_Masonryt(G3_Runtime *rt)
{

  // Pointer to a uniaxial material that will be returned
  // UniaxialMaterial *theMaterial = 0;
  //
  // parse the input line for the material parameters
  //
  int iData[1];
  double dData[21];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Masonryt tag" << endln;
    return 0;
  }
  numData = 21;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Masonryt Material Parameters\n";
    return 0;
  }
  //
  // create a new material
  //
  UniaxialMaterial *mat =
      new Masonryt(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                   dData[5], dData[6], dData[7], dData[8], dData[9], dData[10],
                   dData[11], dData[12], dData[13], dData[14], dData[15],
                   dData[16], dData[17], dData[18], dData[19], int(dData[20]));

  if (mat == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Masonryt\n";
    return 0;
  }
  // return the material
  return mat;
}
