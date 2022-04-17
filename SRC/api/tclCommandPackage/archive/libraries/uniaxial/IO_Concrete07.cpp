

#include <g3_api.h>


#include <SRC/material/uniaxial/Concrete07.h>
void *OPS_Concrete07(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 9) {
    opserr << "WARNING: Insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Concrete07 tag? ";
    opserr << "fpc? epsc0? Ec? fpt? epst0? xcrp? xcrn? r?\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[8];
  numdata = 8;
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double data\n";
    return 0;
  }

  UniaxialMaterial *mat =
      new Concrete07(tag, data[0], data[1], data[2], data[3], data[4], data[5],
                     data[6], data[7]);
  if (mat == 0) {
    opserr << "WARNING: failed to create Concrete07 material\n";
    return 0;
  }

  return mat;
}
