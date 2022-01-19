

#include <g3_api.h>


#include <SRC/material/uniaxial/Trilinwpd.h>
void *OPS_Trilinwpd(G3_Runtime *rt)
{
  // print out some KUDO's
  if (numTrilinwpd == 0) {
    opserr << "Trilineal with pinching unaxial material - Written by GST "
              "UNcuyo Copyright 2017 - Use at your Own Peril\n";
    numTrilinwpd = 1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int iData[2];
  double dData[19];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial trilinwpd tag" << endln;
    return 0;
  }

  numData = 19;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid parameters\n";
    return 0;
  }
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial trilinwpd type" << endln;
    return 0;
  }

  int itype = iData[1];

  //
  // create a new material
  //

  theMaterial =
      new trilinwpd(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                    dData[5], dData[6], dData[7], dData[8], dData[9], dData[10],
                    dData[11], dData[12], dData[13], dData[14], dData[15],
                    dData[16], dData[17], dData[18], iData[1]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type trilinwpd\n";
    return 0;
  }

  // return the material
  return theMaterial;
}
