

#include <g3_api.h>


#include <SRC/material/uniaxial/HardeningMaterial.h>
void *OPS_HardeningMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Hardening tag? E? sigmaY? H_iso? H_kin? "
              "<eta?>"
           << endln;
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING: failed to read tag\n";
    return 0;
  }

  double data[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARING: failed to read data\n";
    return 0;
  }

  double eta = 0.0;
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 0) {
    numdata = 1;
    if (OPS_GetDouble(&numdata, &eta) < 0) {
      opserr << "WARNING: failed to read eta\n";
      return 0;
    }
  }

  UniaxialMaterial *mat =
      new HardeningMaterial(tag, data[0], data[1], data[2], data[3], eta);
  if (mat == 0) {
    opserr << "WARNING: failed to create Hardeningmaterial material\n";
    return 0;
  }

  return mat;
}
