

#include <g3_api.h>


#include <SRC/material/uniaxial/strength/PetrangeliStrengthDegradation.h>
void *OPS_PetrangeliStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Petrangeli "
              "tag? e1? V2? e2?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[3];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Petrangeli" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Petrangeli"
           << endln;
    return 0;
  }

  theDegradation =
      new PetrangeliStrengthDegradation(iData[0], dData[0], dData[1], dData[2]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create PetrangeliStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}
