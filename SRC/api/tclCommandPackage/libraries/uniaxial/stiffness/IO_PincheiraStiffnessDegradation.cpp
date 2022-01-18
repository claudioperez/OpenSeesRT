

#include <g3_api.h>


#include <SRC/material/uniaxial/stiffness/PincheiraStiffnessDegradation.h>
void *OPS_PincheiraStiffnessDegradation(void)
{
  StiffnessDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: stiffnessDegradation Pincheira "
              "tag? alpha? beta? eta? nu?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[4];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for stiffnessDegradation Pincheira" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for stiffnessDegradation Pincheira"
           << endln;
    return 0;
  }

  theDegradation = new PincheiraStiffnessDegradation(
      iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create PincheiraStiffnessDegradation\n";
    return 0;
  }

  return theDegradation;
}
