

#include <g3_api.h>


#include <SRC/analysis/integrator/AlphaOS_TP.h>
void *OPS_AlphaOS_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 1 || argc > 4) {
    opserr << "WARNING - incorrect number of args want AlphaOS_TP $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
    return 0;
  }

  bool updElemDisp = false;
  double dData[3];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 3;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr
        << "WARNING - invalid args want AlphaOS_TP $alpha <-updateElemDisp>\n";
    opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
    return 0;
  }

  if (argc == 2 || argc == 4) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (argc < 3)
    theIntegrator = new AlphaOS_TP(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOS_TP(dData[0], dData[1], dData[2], updElemDisp);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating AlphaOS_TP integrator\n";

  return theIntegrator;
}
