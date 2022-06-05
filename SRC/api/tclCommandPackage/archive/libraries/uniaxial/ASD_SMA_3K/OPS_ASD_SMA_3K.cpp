

#include <g3_api.h>


#include <SRC/material/uniaxial/ASD_SMA_3K.h>
void *OPS_ASD_SMA_3K(G3_Runtime *rt)
{

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING: Insufficient arguments\n";
    opserr
        << "Want: uniaxialMaterial ASD_SMA_3K matTag? k1? k2? k3? sigF? beta?";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[5] = {0, 0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 5) {
    numdata = 5;
  }
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  UniaxialMaterial *mat =
      new ASD_SMA_3K(tag, data[0], data[1], data[2], data[3], data[4]);
  if (mat == 0) {
    opserr << "WARNING: failed to create ASD_SMA_3K material\n";
    return 0;
  }

  return mat;
}
