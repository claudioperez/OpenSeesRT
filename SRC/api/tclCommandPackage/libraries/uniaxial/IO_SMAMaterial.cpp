

#include <g3_api.h>


#include <SRC/material/uniaxial/SMAMaterial.h>
void *OPS_SMAMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial SMA matTag? E? eps_L? sig_AM_s? "
              "sig_AM_f? sig_MA_s? sig_MA_f?"
           << endln;
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING: failed to read tag\n";
    return 0;
  }

  double data[6];
  numdata = 6;
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING: failed to read data\n";
    return 0;
  }

  UniaxialMaterial *mat = new SMAMaterial(tag, data[0], data[1], data[2],
                                          data[3], data[4], data[5]);
  if (mat == 0) {
    opserr << "WARNING: failed to create SMAMaterial\n";
    return 0;
  }

  return mat;
}
