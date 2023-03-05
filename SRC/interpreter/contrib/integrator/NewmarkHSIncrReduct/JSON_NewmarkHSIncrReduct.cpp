

#include <g3_api.h>


#include <SRC/analysis/integrator/NewmarkHSIncrReduct.h>
void *OPS_NewmarkHSIncrReduct(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 3) {
    opserr << "WARNING - incorrect number of args want NewmarkHSIncrReduct "
              "$gamma $beta $reduct\n";
    return 0;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSIncrReduct $gamma $beta "
              "$reduct\n";
    return 0;
  }

  theIntegrator = new NewmarkHSIncrReduct(dData[0], dData[1], dData[2]);

  if (theIntegrator == 0)
    opserr
        << "WARNING - out of memory creating NewmarkHSIncrReduct integrator\n";

  return theIntegrator;
}
