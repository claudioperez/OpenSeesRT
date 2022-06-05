

#include <g3_api.h>


#include <SRC/analysis/integrator/CollocationHSIncrReduct.h>
void *OPS_CollocationHSIncrReduct(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want CollocationHSIncrReduct "
              "$theta $reduct\n";
    opserr
        << "          or CollocationHSIncrReduct $theta $beta $gamma $reduct\n";
    return 0;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSIncrReduct $theta "
              "$reduct\n";
    opserr
        << "          or CollocationHSIncrReduct $theta $beta $gamma $reduct\n";
    return 0;
  }

  if (argc == 2)
    theIntegrator = new CollocationHSIncrReduct(dData[0], dData[1]);
  else
    theIntegrator =
        new CollocationHSIncrReduct(dData[0], dData[1], dData[2], dData[3]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CollocationHSIncrReduct "
              "integrator\n";

  return theIntegrator;
}
