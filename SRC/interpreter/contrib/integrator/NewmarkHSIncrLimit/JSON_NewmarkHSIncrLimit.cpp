

#include <g3_api.h>


#include <SRC/analysis/integrator/NewmarkHSIncrLimit.h>
void *OPS_NewmarkHSIncrLimit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 3 && argc != 5) {
    opserr << "WARNING - incorrect number of args want NewmarkHSIncrLimit "
              "$gamma $beta $limit <-normType $T>\n";
    return 0;
  }

  double dData[3];
  int normType = 2;
  int numData = 3;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSIncrLimit $gamma $beta "
              "$limit <-normType $T>\n";
    return 0;
  }

  if (argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &normType) != 0) {
        opserr << "WARNING - invalid normType want NewmarkHSIncrLimit $gamma "
                  "$beta $limit <-normType $T>\n";
      }
    }
  }

  theIntegrator =
      new NewmarkHSIncrLimit(dData[0], dData[1], dData[2], normType);

  if (theIntegrator == 0)
    opserr
        << "WARNING - out of memory creating NewmarkHSIncrLimit integrator\n";

  return theIntegrator;
}
