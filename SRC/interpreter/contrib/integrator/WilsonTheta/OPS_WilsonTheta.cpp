

#include <g3_api.h>


#include <SRC/analysis/integrator/WilsonTheta.h>
void *OPS_WilsonTheta(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1) {
    opserr << "WARNING - incorrect number of args want WilsonTheta $theta\n";
    return 0;
  }

  double theta;
  if (OPS_GetDouble(&argc, &theta) != 0) {
    opserr << "WARNING - invalid args want WilsonTheta $theta\n";
    return 0;
  }

  theIntegrator = new WilsonTheta(theta);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating WilsonTheta integrator\n";

  return theIntegrator;
}
