

#include <g3_api.h>


#include <SRC/material/uniaxial/strength/ACIStrengthDegradation.h>
void *OPS_ACIStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation ACI tag? Ky? "
              "D1? v2? D2?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[4];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation ACI" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation ACI" << endln;
    return 0;
  }

  theDegradation = new ACIStrengthDegradation(iData[0], dData[0], dData[1],
                                              dData[2], dData[3]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create ACIStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}
