

#include <g3_api.h>


#include <SRC/analysis/integrator/Newmark.h>
void *OPS_Newmark(void) {
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want Newmark $gamma $beta "
              "<-form $typeUnknown>\n";
    return 0;
  }

  int dispFlag = 1;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want Newmark $gamma $beta <-form "
              "$typeUnknown>\n";
    return 0;
  }

  if (argc == 2)
    theIntegrator = new Newmark(dData[0], dData[1]);
  else {
    //    char nextString[10];
    const char *nextString = OPS_GetString();
    //    OPS_GetString(nextString, 10);
    if (strcmp(nextString, "-form") == 0) {
      //      OPS_GetString(nextString, 10);
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd'))
        dispFlag = 1;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        dispFlag = 3;
      else if ((nextString[0] == 'V') || (nextString[0] == 'v'))
        dispFlag = 2;
    }
    theIntegrator = new Newmark(dData[0], dData[1], dispFlag);
  }

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating Newmark integrator\n";

  return theIntegrator;
}
