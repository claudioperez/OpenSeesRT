

#include <g3_api.h>


#include <SRC/analysis/integrator/KRAlphaExplicit_TP.h>
void *OPS_KRAlphaExplicit_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1) {
    opserr << "WARNING - incorrect number of args want KRAlphaExplicit_TP "
              "$rhoInf\n";
    return 0;
  }

  double rhoInf;
  if (OPS_GetDouble(&argc, &rhoInf) != 0) {
    opserr << "WARNING - invalid args want KRAlphaExplicit_TP $rhoInf\n";
    return 0;
  }

  theIntegrator = new KRAlphaExplicit_TP(rhoInf);

  if (theIntegrator == 0)
    opserr
        << "WARNING - out of memory creating KRAlphaExplicit_TP integrator\n";

  return theIntegrator;
}
