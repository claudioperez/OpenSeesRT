

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTHSFixedNumIter_TP.h>
void *OPS_HHTHSFixedNumIter_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 3 && argc != 4 && argc != 6) {
    opserr << "WARNING - incorrect number of args want HHTHSFixedNumIter_TP "
              "$rhoInf <-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return 0;
  }

  double dData[4];
  int polyOrder = 2;
  bool updDomFlag = true;
  int numData;
  if (argc < 4)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSFixedNumIter_TP $rhoInf "
              "<-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return 0;
  }

  if (argc == 3 || argc == 6) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want HHTHSFixedNumIter_TP "
                  "$rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta "
                  "$gamma <-polyOrder $O>\n";
      }
    }
  }

  if (argc < 4)
    theIntegrator = new HHTHSFixedNumIter_TP(dData[0], polyOrder, updDomFlag);
  else
    theIntegrator = new HHTHSFixedNumIter_TP(dData[0], dData[1], dData[2],
                                             dData[3], polyOrder, updDomFlag);

  if (theIntegrator == 0)
    opserr
        << "WARNING - out of memory creating HHTHSFixedNumIter_TP integrator\n";

  return theIntegrator;
}
