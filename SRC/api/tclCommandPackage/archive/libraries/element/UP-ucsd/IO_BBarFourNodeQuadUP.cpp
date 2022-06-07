

#include <g3_api.h>


#include <SRC/element/UP-ucsd/BBarFourNodeQuadUP.h>
void *OPS_BBarFourNodeQuadUP()
{
  if (OPS_GetNDM() != 2 || OPS_GetNDF() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element bbarQuadUP eleTag? iNode? jNode? kNode? lNode? "
              "thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? "
              "pressure? dM? dK?>\n";
    return 0;
  }

  // BBarFourNodeQuadUPId, iNode, jNode, kNode, lNode
  int tags[5];
  int num = 5;
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
    opserr << "Material: " << matTag;
    opserr << "\nBBarFourNodeQuadUP element: " << tags[0] << endln;
    return 0;
  }

  // bk, r, perm1, perm2
  double data[4];
  num = 4;
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: invalid double input\n";
    return 0;
  }

  // b1, b2, p
  double opt[3] = {0, 0, 0};
  num = OPS_GetNumRemainingInputArgs();
  if (num > 3) {
    num = 3;
  }
  if (num > 0) {
    if (OPS_GetDoubleInput(&num, opt) < 0) {
      opserr << "WARNING: invalid double input\n";
      return 0;
    }
  }

  return new BBarFourNodeQuadUP(tags[0], tags[1], tags[2], tags[3], tags[4],
                                *mat, "PlaneStrain", thk, data[0], data[1],
                                data[2], data[3], opt[0], opt[1], opt[2]);
}
