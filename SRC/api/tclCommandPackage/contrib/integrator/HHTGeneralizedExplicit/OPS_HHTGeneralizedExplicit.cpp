

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTGeneralizedExplicit.h>
void *OPS_HHTGeneralizedExplicit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2 || argc > 5) {
    opserr << "WARNING - incorrect number of args want HHTGeneralizedExplicit "
              "$rhoB $alphaF <-updateElemDisp>\n";
    opserr << "          or HHTGeneralizedExplicit $alphaI $alphaF $beta "
              "$gamma <-updateElemDisp>\n";
    return 0;
  }

  bool updElemDisp = false;
  double dData[4];
  int numData;
  if (argc < 4)
    numData = 2;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralizedExplicit $rhoB "
              "$alphaF <-updateElemDisp>\n";
    opserr << "          or HHTGeneralizedExplicit $alphaI $alphaF $beta "
              "$gamma <-updateElemDisp>\n";
    return 0;
  }

  if (argc == 3 || argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (argc < 4)
    theIntegrator = new HHTGeneralizedExplicit(dData[0], dData[1], updElemDisp);
  else
    theIntegrator = new HHTGeneralizedExplicit(dData[0], dData[1], dData[2],
                                               dData[3], updElemDisp);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTGeneralizedExplicit "
              "integrator\n";

  return theIntegrator;
}
