

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTExplicit.h>
void *OPS_HHTExplicit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 1 || argc > 3) {
    opserr << "WARNING - incorrect number of args want HHTExplicit $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
    return 0;
  }

  bool updElemDisp = false;
  double dData[2];
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr
        << "WARNING - invalid args want HHTExplicit $alpha <-updateElemDisp>\n";
    opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
    return 0;
  }

  if (numData + 1 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (numData == 1)
    theIntegrator = new HHTExplicit(dData[0], updElemDisp);
  else if (numData == 2)
    theIntegrator = new HHTExplicit(dData[0], dData[1], updElemDisp);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTExplicit integrator\n";

  return theIntegrator;
}
