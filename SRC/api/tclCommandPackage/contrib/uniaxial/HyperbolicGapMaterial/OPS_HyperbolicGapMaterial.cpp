

#include <g3_api.h>


#include <SRC/material/uniaxial/HyperbolicGapMaterial.h>
void *OPS_HyperbolicGapMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 6) {
    opserr << "WARNING: Insufficient arguments\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    return 0;
  }

  double data[5];
  numdata = 5;
  if (OPS_GetDoubleInput(&numdata, data)) {
    return 0;
  }

  UniaxialMaterial *mat = new HyperbolicGapMaterial(tag, data[0], data[1],
                                                    data[2], data[3], data[4]);
  if (mat == 0) {
    opserr << "WARNING: failed to create Hyperbolicgapmaterial material\n";
    return 0;
  }

  return mat;
}
