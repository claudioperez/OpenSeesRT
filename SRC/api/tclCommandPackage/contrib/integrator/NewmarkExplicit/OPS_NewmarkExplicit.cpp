

#include <g3_api.h>


#include <SRC/analysis/integrator/NewmarkExplicit.h>
void *OPS_NewmarkExplicit(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1) {
    opserr
        << "WARNING - incorrect number of args want NewmarkExplicit $gamma\n";
    return 0;
  }

  double gamma;
  if (OPS_GetDouble(&argc, &gamma) != 0) {
    opserr << "WARNING - invalid args want NewmarkExplicit $gamma\n";
    return 0;
  }

  theIntegrator = new NewmarkExplicit(gamma);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating NewmarkExplicit integrator\n";

  return theIntegrator;
}
