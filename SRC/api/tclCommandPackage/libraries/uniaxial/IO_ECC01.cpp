

#include <g3_api.h>


#include <SRC/material/uniaxial/ECC01.h>
void *OPS_ECC01(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial ECC01 TAG? SIGT0? EPST0? SIGT1? EPST1? "
              "EPST2? SIGC0? EPSC0? EPSC1? ";
    opserr << "ALPHAT1? ALPHAT2? ALPHAC? ALPHACU? BETAT? BETAC\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    return 0;
  }

  double data[14];
  numdata = 14;
  if (OPS_GetDoubleInput(&numdata, data)) {
    return 0;
  }

  UniaxialMaterial *mat = new ECC01(
      tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6],
      data[7], data[8], data[9], data[10], data[11], data[12], data[13]);
  if (mat == 0) {
    opserr << "WARNING: failed to create ECC01 material\n";
    return 0;
  }

  return mat;
}
