

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTExplicit_TP.h>
void *OPS_HHTExplicit_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 1 || argc > 2) {
    opserr << "WARNING - incorrect number of args want HHTExplicit_TP $alpha\n";
    opserr << "          or HHTExplicit_TP $alpha $gamma\n";
    return 0;
  }

  double dData[2];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTExplicit_TP $alpha\n";
    opserr << "          or HHTExplicit_TP $alpha $gamma\n";
    return 0;
  }

  if (argc == 1)
    theIntegrator = new HHTExplicit_TP(dData[0]);
  else if (argc == 2)
    theIntegrator = new HHTExplicit_TP(dData[0], dData[1]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTExplicit_TP integrator\n";

  return theIntegrator;
}
