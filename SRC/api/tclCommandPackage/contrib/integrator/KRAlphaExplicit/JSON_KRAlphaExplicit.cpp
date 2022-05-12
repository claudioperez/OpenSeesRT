

#include <g3_api.h>


#include <SRC/analysis/integrator/KRAlphaExplicit.h>
void *OPS_KRAlphaExplicit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 2) {
    opserr << "WARNING - incorrect number of args want KRAlphaExplicit $rhoInf "
              "<-updateElemDisp>\n";
    return 0;
  }

  bool updElemDisp = false;
  double rhoInf;
  int numData = 1;
  if (OPS_GetDouble(&numData, &rhoInf) != 0) {
    opserr << "WARNING - invalid args want KRAlphaExplicit $rhoInf "
              "<-updateElemDisp>\n";
    return 0;
  }

  if (argc == 2) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  theIntegrator = new KRAlphaExplicit(rhoInf, updElemDisp);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating KRAlphaExplicit integrator\n";

  return theIntegrator;
}
