

#include <g3_api.h>


#include <SRC/analysis/integrator/CentralDifferenceNoDamping.h>
void *OPS_CentralDifferenceNoDamping(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  theIntegrator = new CentralDifferenceNoDamping();

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CentralDifferenceNoDamping "
              "integrator\n";

  return theIntegrator;
}
