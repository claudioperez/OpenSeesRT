

#include <g3_api.h>


#include <SRC/analysis/integrator/CollocationHSIncrLimit.h>
void *OPS_CollocationHSIncrLimit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4 && argc != 6) {
    opserr << "WARNING - incorrect number of args want CollocationHSIncrLimit "
              "$theta $limit <-normType $T>\n";
    opserr << "          or CollocationHSIncrLimit $theta $beta $gamma $limit "
              "<-normType $T>\n";
    return 0;
  }

  double dData[4];
  int normType = 2;
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSIncrLimit $theta "
              "$limit <-normType $T>\n";
    opserr << "          or CollocationHSIncrLimit $theta $beta $gamma $limit "
              "<-normType $T>\n";
    return 0;
  }

  if (numData + 2 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      int numData2 = 1;
      if (OPS_GetInt(&numData2, &normType) != 0) {
        opserr << "WARNING - invalid normType want CollocationHSIncrLimit "
                  "$theta $limit <-normType $T>\n";
        opserr << "          or CollocationHSIncrLimit $theta $beta $gamma "
                  "$limit <-normType $T>\n";
      }
    }
  }

  if (numData == 2)
    theIntegrator = new CollocationHSIncrLimit(dData[0], dData[1], normType);
  else if (numData == 4)
    theIntegrator = new CollocationHSIncrLimit(dData[0], dData[1], dData[2],
                                               dData[3], normType);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CollocationHSIncrLimit "
              "integrator\n";

  return theIntegrator;
}
