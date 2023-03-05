

#include <g3_api.h>


#include <SRC/analysis/integrator/AlphaOSGeneralized_TP.h>
void *OPS_AlphaOSGeneralized_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 2 && argc != 4 && argc != 5) {
    opserr << "WARNING - incorrect number of args want AlphaOSGeneralized_TP "
              "$rhoInf <-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized_TP $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return 0;
  }

  bool updElemDisp = false;
  double dData[4];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want AlphaOSGeneralized_TP $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized_TP $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return 0;
  }

  if (argc == 2 || argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (argc < 3)
    theIntegrator = new AlphaOSGeneralized_TP(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOSGeneralized_TP(dData[0], dData[1], dData[2],
                                              dData[3], updElemDisp);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating AlphaOSGeneralized_TP "
              "integrator\n";

  return theIntegrator;
}
