

#include <g3_api.h>


#include <SRC/material/uniaxial/SelfCenteringMaterial.h>
void *OPS_SelfCenteringMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING: Insufficient arguments\n";
    opserr << "Want: uniaxialMaterial SelfCentering tag? k1? k2? ";
    opserr << "ActF? beta? <SlipDef? BearDef? rBear?>" << endln;
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[7] = {0, 0, 0, 0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 7) {
    numdata = 7;
  }
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  UniaxialMaterial *mat = new SelfCenteringMaterial(
      tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6]);
  if (mat == 0) {
    opserr << "WARNING: failed to create Selfcenteringmaterial material\n";
    return 0;
  }

  return mat;
}
