

#include <g3_api.h>


#include <SRC/analysis/integrator/PFEMIntegrator.h>
void *OPS_PFEMIntegrator(void) {
  TransientIntegrator *theIntegrator = 0;

  int dispFlag = 2;
  int init = 2;
  double dData[2] = {-1.0, -1.0};
  int numData = 2;
  if (OPS_GetNumRemainingInputArgs() > 1) {
    if (OPS_GetDouble(&numData, dData) < 0) {
      OPS_ResetCurrentInputArg(-2);
    }
  }

  if (OPS_GetNumRemainingInputArgs() > 1) {
    const char *nextString = OPS_GetString();
    if (strcmp(nextString, "-form") == 0) {
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd')) {
        dispFlag = 1;
        init = 1;
      } else if ((nextString[0] == 'A') || (nextString[0] == 'a')) {
        dispFlag = 3;
        init = 3;
      } else if ((nextString[0] == 'V') || (nextString[0] == 'v')) {
        dispFlag = 2;
        init = 2;
      }
    } else {
      opserr << "WARNING: first option must be -form\n";
      return 0;
    }
    if (OPS_GetNumRemainingInputArgs() > 1) {
      nextString = OPS_GetString();
      if (strcmp(nextString, "-init") == 0) {
        nextString = OPS_GetString();
        if ((nextString[0] == 'D') || (nextString[0] == 'd')) {
          init = 1;
        } else if ((nextString[0] == 'A') || (nextString[0] == 'a')) {
          init = 3;
        } else if ((nextString[0] == 'V') || (nextString[0] == 'v')) {
          init = 2;
        }
      }
    } else {
      opserr << "WARNING: second option must be -init\n";
      return 0;
    }
  }

  theIntegrator = new PFEMIntegrator(dData[0], dData[1], dispFlag, init);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating Newmark integrator\n";

  return theIntegrator;
}
