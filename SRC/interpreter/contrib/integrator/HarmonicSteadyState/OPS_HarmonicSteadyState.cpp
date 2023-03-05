

#include <g3_api.h>


#include <SRC/analysis/integrator/HarmonicSteadyState.h>
void *OPS_HarmonicSteadyState() {
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "insufficient arguments\n";
    return 0;
  }

  double lambda;
  int numData = 1;
  if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
    opserr << "WARNING failed to read double lambda\n";
    return 0;
  }

  double period = 0;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &period) < 0) {
    opserr << "WARNING failed to read double period\n";
    return 0;
  }

  int numIter = 1;
  double mLambda[2] = {lambda, lambda};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING failed to read int numIter\n";
      return 0;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
      opserr << "WARNING failed to read double min and max\n";
      return 0;
    }
  }

  return new HarmonicSteadyState(lambda, period, numIter, mLambda[0],
                                 mLambda[1]);
}
