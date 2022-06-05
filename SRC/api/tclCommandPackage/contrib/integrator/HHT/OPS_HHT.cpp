#include <g3_api.h>
#include <SRC/analysis/integrator/HHT.h>

void *OPS_HHT(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 3) {
    opserr << "WARNING - incorrect number of args want HHT $alpha <$gamma "
              "$beta>\n";
    return 0;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHT $alpha <$gamma $beta>\n";
    return 0;
  }

  if (argc == 1)
    theIntegrator = new HHT(dData[0]);
  else
    theIntegrator = new HHT(dData[0], dData[1], dData[2]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHT integrator\n";

  return theIntegrator;
}
