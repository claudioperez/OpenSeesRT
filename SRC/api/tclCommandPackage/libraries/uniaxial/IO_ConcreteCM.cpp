

#include <g3_api.h>


#include <SRC/material/uniaxial/ConcreteCM.h>
void *OPS_ConcreteCM(void)
{

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  // Parse the script for material parameters
  if (numArgs != 10 && numArgs != 11 && numArgs != 12) {
    opserr << "Incorrect # args Want: uniaxialMaterial ConcreteCM tag? fpcc? "
              "epcc? Ec? rc? xcrn? ft? et? rt? xcrp? <-GapClose gap?>"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[9];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial ConcreteCM ConcreteCM"
           << endln;
    return 0;
  }

  numData = 9;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxialMaterial ConcreteCM ConcreteCM "
           << iData[0] << endln;
    return 0;
  }

  if (numArgs == 10) {

    theMaterial =
        new ConcreteCM(iData[0], dData[0], dData[1], dData[2], dData[3],
                       dData[4], dData[5], dData[6], dData[7], dData[8]);

  } else if (numArgs == 11) {

    numData = 1;
    int mon;
    if (OPS_GetIntInput(&numData, &mon) != 0) {
      opserr
          << "Invalid $mon parameter for uniaxialMaterial ConcreteCM with tag  "
          << iData[0] << endln;
      return 0;
    }

    if (mon != 0 && mon != 1) {
      opserr
          << "Invalid $mon parameter for uniaxialMaterial ConcreteCM with tag  "
          << iData[0] << endln;
      return 0;
    }

    theMaterial =
        new ConcreteCM(iData[0], dData[0], dData[1], dData[2], dData[3],
                       dData[4], dData[5], dData[6], dData[7], dData[8], mon);

  } else {

    int gap;
    numData = 1;

    const char *str = OPS_GetString();
    // OPS_GetStringCopy(&str);
    if (strcmp(str, "-GapClose") == 0) {
      if (OPS_GetIntInput(&numData, &gap) != 0) {
        opserr << "Invalid $gap parameter for uniaxialMaterial ConcreteCM with "
                  "tag  "
               << iData[0] << endln;
        return 0;
      }
    } else {
      opserr << "Invalid input parameter for uniaxialMaterial ConcreteCM with "
                "tag  "
             << iData[0] << ", want: -GapClose" << endln;
      return 0;
    }

    // delete [] str;

    if (gap != 0 && gap != 1) {
      opserr
          << "Invalid $gap parameter for uniaxialMaterial ConcreteCM with tag  "
          << iData[0] << endln;
      return 0;
    }

    int check = 0; // dummy variable
    theMaterial = new ConcreteCM(iData[0], dData[0], dData[1], dData[2],
                                 dData[3], dData[4], dData[5], dData[6],
                                 dData[7], dData[8], gap, check);
  }

  return theMaterial;
}
