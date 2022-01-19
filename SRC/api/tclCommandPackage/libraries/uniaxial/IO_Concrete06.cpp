

#include <g3_api.h>


#include <SRC/material/uniaxial/Concrete06.h>
void *OPS_Concrete06(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Concrete06 ";
    opserr << "tag? fc? eo? r? k? alphaC? fcr? ecr? b? alphaT?\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[9];
  numdata = 9;
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double data\n";
    return 0;
  }

  UniaxialMaterial *mat =
      new Concrete06(tag, data[0], data[1], data[2], data[3], data[4], data[5],
                     data[6], data[7], data[8]);
  if (mat == 0) {
    opserr << "WARNING: failed to create Concrete06 material\n";
    return 0;
  }

  return mat;
}
