

#include <g3_api.h>


#include <SRC/material/uniaxial/unloading/EnergyUnloadingRule.h>
void *OPS_EnergyUnloadingRule(void)
{
  UnloadingRule *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: unloadingRule Energy tag? Et? c?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[2];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for unloadingRule Energy" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for unloadingRule Energy" << endln;
    return 0;
  }

  theDegradation = new EnergyUnloadingRule(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create EnergyUnloadingRule\n";
    return 0;
  }

  return theDegradation;
}
