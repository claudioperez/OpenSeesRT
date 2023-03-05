

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTGeneralized.h>
void *OPS_HHTGeneralized(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 4) {
    opserr
        << "WARNING - incorrect number of args want HHTGeneralized $rhoInf\n";
    opserr << "          or HHTGeneralized $alphaI $alphaF $beta $gamma\n";
    return 0;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralized $rhoInf\n";
    opserr << "          or HHTGeneralized $alphaI $alphaF $beta $gamma\n";
    return 0;
  }

  if (argc == 1)
    theIntegrator = new HHTGeneralized(dData[0]);
  else
    theIntegrator = new HHTGeneralized(dData[0], dData[1], dData[2], dData[3]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTGeneralized integrator\n";

  return theIntegrator;
}
