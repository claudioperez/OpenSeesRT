

#include <g3_api.h>


#include <SRC/material/uniaxial/BoucWenMaterial.h>
void *OPS_BoucWenMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 10) {
    opserr << "WARNING: Insufficient arguments\n";
    opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?"
           << endln << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0e-8};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 10) {
    numdata = 10;
  }
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  int maxNumIter = 20;
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 0) {
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &maxNumIter) < 0) {
      opserr << "WARNING invalid int inputs\n";
      return 0;
    }
  }

  UniaxialMaterial *mat = new BoucWenMaterial(
      tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6],
      data[7], data[8], data[9], maxNumIter);
  if (mat == 0) {
    opserr << "WARNING: failed to create Boucwenmaterial material\n";
    return 0;
  }

  return mat;
}
