

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTHSIncrReduct.h>
void *OPS_HHTHSIncrReduct(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 5) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrReduct $rhoInf "
              "$reduct\n";
    opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return 0;
  }

  double dData[5];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSIncrReduct $rhoInf $reduct\n";
    opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return 0;
  }

  if (argc == 2)
    theIntegrator = new HHTHSIncrReduct(dData[0], dData[1]);
  else
    theIntegrator =
        new HHTHSIncrReduct(dData[0], dData[1], dData[2], dData[3], dData[4]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTHSIncrReduct integrator\n";

  return theIntegrator;
}
