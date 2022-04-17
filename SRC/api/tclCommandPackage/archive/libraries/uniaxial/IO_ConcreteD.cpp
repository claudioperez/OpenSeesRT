

#include <g3_api.h>


#include <SRC/material/uniaxial/ConcreteD.h>
void *OPS_ConcreteD(void)
{
  // print out some KUDO's
  if (numConcreteD == 0) {
    // opserr << "ConcreteD unaxial material - Written by Zenyong Wan Tongji
    // University 2014.01\n";
    numConcreteD = 1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int iData[1];
  double dData[9];

  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ConcreteD tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 7 && numData != 9) {
    opserr << "Invalid #args, want: uniaxialMaterial ConcreteD " << iData[0]
           << "(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)"
           << endln;
    return 0;
  }

  if (numData == 7) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args: uniaxialMaterial ConcreteD " << iData[0]
             << "(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)"
             << endln;
      return 0;
    }
    // Parsing was successful, allocate the material
    theMaterial = new ConcreteD(iData[0], dData[0], dData[1], dData[2],
                                dData[3], dData[4], dData[5], dData[6]);

  } else if (numData == 9) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args: uniaxialMaterial ConcreteD " << iData[0]
             << "(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)"
             << endln;
      return 0;
    }
    // Parsing was successful, allocate the material
    theMaterial =
        new ConcreteD(iData[0], dData[0], dData[1], dData[2], dData[3],
                      dData[4], dData[5], dData[6], dData[7], dData[8]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteD\n";
    return 0;
  }

  // return the material
  return theMaterial;
}
