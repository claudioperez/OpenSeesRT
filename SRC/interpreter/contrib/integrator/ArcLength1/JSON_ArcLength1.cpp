

#include <g3_api.h>


#include <SRC/analysis/integrator/ArcLength1.h>
void *OPS_ArcLength1() {
  double arcLength;
  double alpha;
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING integrator ArcLength arcLength alpha \n";
    return 0;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
    opserr << "WARNING integrator ArcLength failed to read arc length\n";
    return 0;
  }
  if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
    opserr << "WARNING integrator ArcLength failed to read alpha\n";
    return 0;
  }
  return new ArcLength1(arcLength, alpha);
}
