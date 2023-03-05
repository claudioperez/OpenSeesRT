

#include <g3_api.h>


#include <SRC/analysis/integrator/Collocation.h>
void *OPS_Collocation(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 1 && argc != 3) {
    opserr << "WARNING - incorrect number of args want Collocation $theta\n";
    opserr << "          or Collocation $theta $beta $gamma\n";
    return 0;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want Collocation $theta\n";
    opserr << "          or Collocation $theta $beta $gamma\n";
    return 0;
  }

  if (argc == 1)
    theIntegrator = new Collocation(dData[0]);
  else
    theIntegrator = new Collocation(dData[0], dData[1], dData[2]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating Collocation integrator\n";

  return theIntegrator;
}
