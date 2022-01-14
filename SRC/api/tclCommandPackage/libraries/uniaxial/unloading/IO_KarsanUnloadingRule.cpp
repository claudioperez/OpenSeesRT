

#include <g3_api.h>


#include <SRC/material/uniaxial/unloading/KarsanUnloadingRule.h>
void *OPS_KarsanUnloadingRule(void)
{
  UnloadingRule *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr
        << "Invalid number of args, want: unloadingRule Karsan tag? epsc? epsu?"
        << endln;
    return 0;
  }

  int iData[1];
  double dData[2];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for unloadingRule Karsan" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for unloadingRule Karsan" << endln;
    return 0;
  }

  theDegradation = new KarsanUnloadingRule(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create KarsanUnloadingRule\n";
    return 0;
  }

  return theDegradation;
}
