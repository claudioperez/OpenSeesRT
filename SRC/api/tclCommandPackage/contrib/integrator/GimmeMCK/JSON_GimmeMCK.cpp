

#include <g3_api.h>


#include <SRC/analysis/integrator/GimmeMCK.h>
void *OPS_GimmeMCK(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3) {
    opserr
        << "WARNING - incorrect number of args want GimmeMCK $m $c $k <$ki>\n";
    return 0;
  }

  int numdata = 3;
  double ddata[3];
  if (OPS_GetDouble(&numdata, ddata) != 0) {
    opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
    return 0;
  }
  numdata = 1;
  double ki = 0.0;
  if (argc > 3) {
    if (OPS_GetDouble(&numdata, &ki) != 0) {
      opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
      return 0;
    }
  }

  theIntegrator = new GimmeMCK(ddata[0], ddata[1], ddata[2], ki);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating GimmeMCK integrator\n";

  return theIntegrator;
}
