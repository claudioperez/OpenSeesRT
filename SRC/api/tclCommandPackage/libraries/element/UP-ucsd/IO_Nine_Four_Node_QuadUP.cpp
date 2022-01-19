

#include <g3_api.h>


#include <SRC/element/UP-ucsd/Nine_Four_Node_QuadUP.h>
void *OPS_NineFourNodeQuadUP()
{
  if (OPS_GetNDM() != 2) {
    opserr << "WARNING -- model dimensions not compatible with 9-4-NodeQuadUP "
              "element\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 16) {
    opserr << "WARNING insufficient arguments\n";
    opserr
        << "Want: element FourNodeQuadUP eleTag? Node1? ... Node9? thk? type? "
           "matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return 0;
  }

  // NineFourNodeQuadUPId, Node[9]
  int tags[10];
  int num = 10;
  if (OPS_GetIntInput(&num, tags) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }

  double thk;
  num = 1;
  if (OPS_GetDoubleInput(&num, &thk) < 0) {
    opserr << "WARNING: invalid double input\n";
    return 0;
  }

  int matTag;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }
  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << matTag;
    opserr << "\nQuad element: " << tags[0] << endln;
  }

  // bk, r, perm1, perm2
  double data[4];
  num = 4;
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: invalid double input\n";
    return 0;
  }

  // b1, b2
  double opt[2] = {0, 0};
  num = OPS_GetNumRemainingInputArgs();
  if (num > 2) {
    num = 2;
  }
  if (num > 0) {
    if (OPS_GetDoubleInput(&num, opt) < 0) {
      opserr << "WARNING: invalid double input\n";
      return 0;
    }
  }

  return new NineFourNodeQuadUP(tags[0], tags[1], tags[2], tags[3], tags[4],
                                tags[5], tags[6], tags[7], tags[8], tags[9],
                                *mat, "PlaneStrain", thk, data[0], data[1],
                                data[2], data[3], opt[0], opt[1]);
}
