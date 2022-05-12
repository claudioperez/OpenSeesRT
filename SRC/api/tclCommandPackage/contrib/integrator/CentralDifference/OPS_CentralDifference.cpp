

#include <g3_api.h>


#include <SRC/analysis/integrator/CentralDifference.h>
void *OPS_CentralDifference(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  theIntegrator = new CentralDifference();

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CentralDifference integrator\n";

  return theIntegrator;
}
