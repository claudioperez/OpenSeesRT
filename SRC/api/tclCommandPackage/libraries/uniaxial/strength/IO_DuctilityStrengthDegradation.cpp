

#include <g3_api.h>


#include <SRC/material/uniaxial/strength/DuctilityStrengthDegradation.h>
void *OPS_DuctilityStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Ductility "
              "tag? alpha? beta?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[2];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Ductility" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Ductility" << endln;
    return 0;
  }

  theDegradation =
      new DuctilityStrengthDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create DuctilityStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}
