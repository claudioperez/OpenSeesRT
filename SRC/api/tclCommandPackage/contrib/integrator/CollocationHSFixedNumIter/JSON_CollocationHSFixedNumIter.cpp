

#include <g3_api.h>


#include <SRC/analysis/integrator/CollocationHSFixedNumIter.h>
void *OPS_CollocationHSFixedNumIter(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 3 && argc != 5) {
    opserr << "WARNING - incorrect number of args want "
              "CollocationHSFixedNumIter $theta <-polyOrder $O>\n";
    opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma "
              "<-polyOrder $O>\n";
    return 0;
  }

  double dData[3];
  int polyOrder = 2;
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSFixedNumIter $theta "
              "<-polyOrder $O>\n";
    opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma "
              "<-polyOrder $O>\n";
    return 0;
  }

  if (numData + 2 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      int numData2 = 1;
      if (OPS_GetInt(&numData2, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want CollocationHSFixedNumIter "
                  "$rhoInf <-polyOrder $O>\n";
        opserr << "          or CollocationHSFixedNumIter $alphaI $alphaF "
                  "$beta $gamma <-polyOrder $O>\n";
      }
    }
  }

  if (numData == 1)
    theIntegrator = new CollocationHSFixedNumIter(dData[0], polyOrder);
  else if (numData == 3)
    theIntegrator =
        new CollocationHSFixedNumIter(dData[0], dData[1], dData[2], polyOrder);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CollocationHSFixedNumIter "
              "integrator\n";

  return theIntegrator;
}
