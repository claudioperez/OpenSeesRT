

#include <g3_api.h>


#include <SRC/material/uniaxial/Trilinwp2.h>
void *OPS_Trilinwp2(G3_Runtime *rt)
{
  // print out some KUDO's
  // if (numTrilinwp == 0) {
  //  opserr << "Trilineal with pinching unaxial material - Written by GST
  //  UNcuyo Copyright 2017 - Use at your Own Peril\n"; numTrilinwp =1;
  //}

  // Pointer to a uniaxial material that will be returned
  // UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int iData[2];
  double dData[15];
  int numData;
  int numDatatot;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Trilinwp2 tag" << endln;
    return 0;
  }
  numDatatot = numData;

  numData = 15;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid parameters\n";
    return 0;
  }
  numDatatot = numDatatot + numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Trilinwp2 type" << endln;
    return 0;
  }
  numDatatot = numDatatot + numData;
  if (numDatatot != 17) {
    opserr
        << "Invalid Args want: uniaxialMaterial Trilinwp2 tag? Fcrp? dcrp? "
           "Fyp? dyp? Fup? dup? px? py? d1? d2? beta? Pt? Pb? Pc? Mb? itype?  ";
    return 0;
  }

  int itype = iData[1];

  //
  // create a new material
  //

  UniaxialMaterial *mat =
      new Trilinwp2(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                    dData[5], dData[6], dData[7], dData[8], dData[9], dData[10],
                    dData[11], dData[12], dData[13], dData[14], iData[1]);

  if (mat == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Trilinwp2\n";
    return 0;
  }

  // return the material
  return mat;
}
