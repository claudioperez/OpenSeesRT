

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTHSIncrLimit_TP.h>
void *OPS_HHTHSIncrLimit_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4 && argc != 5 && argc != 7) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrLimit_TP "
              "$rhoInf $limit <-normType $T>\n";
    opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
              "$limit <-normType $T>\n";
    return 0;
  }

  double dData[5];
  int normType = 2;
  int numData;
  if (argc < 5)
    numData = 2;
  else
    numData = 5;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSIncrLimit_TP $rhoInf $limit "
              "<-normType $T>\n";
    opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
              "$limit <-normType $T>\n";
    return 0;
  }

  if (argc == 4 || argc == 7) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &normType) != 0) {
        opserr << "WARNING - invalid normType want HHTHSIncrLimit_TP $rhoInf "
                  "$limit <-normType $T>\n";
        opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
                  "$limit <-normType $T>\n";
      }
    }
  }

  if (argc < 5)
    theIntegrator = new HHTHSIncrLimit_TP(dData[0], dData[1], normType);
  else
    theIntegrator = new HHTHSIncrLimit_TP(dData[0], dData[1], dData[2],
                                          dData[3], dData[4], normType);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTHSIncrLimit_TP integrator\n";

  return theIntegrator;
}
