

#include <g3_api.h>


#include <SRC/analysis/integrator/NewmarkHSFixedNumIter.h>
void *OPS_NewmarkHSFixedNumIter(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want NewmarkHSFixedNumIter "
              "$gamma $beta <-polyOrder $O>\n";
    return 0;
  }

  double dData[2];
  int polyOrder = 2;
  bool updDomFlag = true;
  int numData = 2;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSFixedNumIter $gamma $beta "
              "<-polyOrder $O>\n";
    return 0;
  }

  if (argc == 4) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want NewmarkHSFixedNumIter "
                  "$gamma $beta <-polyOrder $O>\n";
      }
    }
  }

  theIntegrator =
      new NewmarkHSFixedNumIter(dData[0], dData[1], polyOrder, updDomFlag);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating NewmarkHSFixedNumIter "
              "integrator\n";

  return theIntegrator;
}
