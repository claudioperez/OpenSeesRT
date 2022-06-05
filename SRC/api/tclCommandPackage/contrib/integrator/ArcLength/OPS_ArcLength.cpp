

#include <g3_api.h>


#include <SRC/analysis/integrator/ArcLength.h>
void *OPS_ArcLength() {
  double arcLength;
  double alpha;
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING integrator ArcLength arcLength alpha \n";
    return 0;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
    opserr << "WARNING integrator ArcLength failed to read arc lenght\n";
    return 0;
  }
  if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
    opserr << "WARNING integrator ArcLength failed to read alpha\n";
    return 0;
  }
  return new ArcLength(arcLength, alpha);
}
