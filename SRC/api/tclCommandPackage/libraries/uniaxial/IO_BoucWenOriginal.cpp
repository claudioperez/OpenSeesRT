

#include <g3_api.h>


#include <SRC/material/uniaxial/BoucWenOriginal.h>
void *OPS_BoucWenOriginal(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING: Insufficient arguments\n";
    opserr << "Want: uniaxialMaterial BoucWenOriginal tag E fy alphaL" << endln;
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[9] = {0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.5, 0.5, 1.0E-8};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 9) {
    numdata = 9;
  }
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  int maxIter = 25;
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 0) {
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
      opserr << "WARNING invalid int inputs\n";
      return 0;
    }
  }

  UniaxialMaterial *mat =
      new BoucWenOriginal(tag, data[0], data[1], data[2], data[3], data[4],
                          data[5], data[6], data[7], data[8], maxIter);
  if (mat == 0) {
    opserr << "WARNING: failed to create BoucWenOriginal material\n";
    return 0;
  }

  return mat;
}
