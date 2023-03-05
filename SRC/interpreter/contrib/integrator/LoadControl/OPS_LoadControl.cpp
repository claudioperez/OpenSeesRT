

#include <g3_api.h>


#include <SRC/analysis/integrator/LoadControl.h>
void *OPS_LoadControlIntegrator() {
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "LoadControl - insufficient arguments\n";
    return 0;
  }

  double lambda;
  int numData = 1;
  if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
    opserr << "WARNING LoadControl - failed to read double lambda\n";
    return 0;
  }

  int numIter = 1;
  double mLambda[2] = {lambda, lambda};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING LoadControl - failed to read int numIter\n";
      return 0;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
      opserr << "WARNING LoadControl - failed to read double min and max\n";
      return 0;
    }
  }

  return new LoadControl(lambda, numIter, mLambda[0], mLambda[1]);
}
