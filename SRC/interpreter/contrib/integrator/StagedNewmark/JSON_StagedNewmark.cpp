

#include <g3_api.h>


#include <SRC/analysis/integrator/StagedNewmark.h>
void *OPS_StagedNewmark(void) {
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want StagedNewmark $gamma "
              "$beta <-form $typeUnknown>\n";
    return 0;
  }

  bool dispFlag = true;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want StagedNewmark $gamma $beta <-form "
              "$typeUnknown>\n";
    return 0;
  }

  if (argc == 2)
    theIntegrator = new StagedNewmark(dData[0], dData[1]);
  else {
    //    char nextString[10];
    const char *nextString = OPS_GetString();
    //    OPS_GetString(nextString, 10);
    if (strcmp(nextString, "-form") == 0) {
      //      OPS_GetString(nextString, 10);
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd'))
        dispFlag = true;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        dispFlag = false;
    }
    theIntegrator = new StagedNewmark(dData[0], dData[1], dispFlag);
  }

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating StagedNewmark integrator\n";

  return theIntegrator;
}
