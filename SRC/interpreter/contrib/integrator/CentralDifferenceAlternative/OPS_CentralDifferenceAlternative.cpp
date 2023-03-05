

#include <g3_api.h>


#include <SRC/analysis/integrator/CentralDifferenceAlternative.h>
void *OPS_CentralDifferenceAlternative(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  theIntegrator = new CentralDifferenceAlternative();

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CentralDifferenceAlternative "
              "integrator\n";

  return theIntegrator;
}
