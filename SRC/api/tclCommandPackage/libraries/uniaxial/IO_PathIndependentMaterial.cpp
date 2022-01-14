

#include <g3_api.h>


#include <SRC/material/uniaxial/PathIndependentMaterial.h>
void *OPS_PathIndependentMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
    return 0;
  }

  int tag[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, tag) < 0) {
    return 0;
  }

  UniaxialMaterial *theMat = OPS_getUniaxialMaterial(tag[1]);
  if (theMat == 0) {
    opserr << "WARNING material does not exist\n";
    opserr << "material: " << tag[1];
    opserr << "\nuniaxialMaterial PathIndependent: " << tag[0] << endln;
    return 0;
  }

  UniaxialMaterial *mat = new PathIndependentMaterial(tag[0], *theMat);
  if (mat == 0) {
    opserr << "WARNING: failed to create PathIndependentmaterial material\n";
    return 0;
  }

  return mat;
}
