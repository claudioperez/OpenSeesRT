

#include <g3_api.h>


#include <SRC/material/uniaxial/ENTMaterial.h>
void *OPS_ENTMaterial(G3_Runtime *rt)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING: invalid #args: ENT matTag E\n";
    return 0;
  }

  int tag;
  int num = 1;
  if (OPS_GetIntInput(&num, &tag) < 0)
    return 0;

  double E;
  if (OPS_GetDoubleInput(&num, &E) < 0)
    return 0;

  UniaxialMaterial *mat = new ENTMaterial(tag, E);
  if (mat == 0)
    return 0;

  // if(OPS_addUniaxialMaterial(mat) == false) {
  // 	opserr<<"WARNING: failed to add ENT material\n";
  // 	delete mat;
  // 	return 0;
  // }
  return mat;
}
