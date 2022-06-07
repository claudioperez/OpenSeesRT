

#include <g3_api.h>


#include <SRC/analysis/integrator/ExplicitDifference.h>
void *OPS_ExplicitDifference(void) {
  TransientIntegrator *theIntegrator = 0;
  theIntegrator = new ExplicitDifference();

  if (theIntegrator == 0)
    opserr
        << "WARNING - out of memory creating ExplicitDifference integrator\n";

  return theIntegrator;
}
