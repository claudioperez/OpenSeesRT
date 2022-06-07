

#include <g3_api.h>


#include <SRC/element/UP-ucsd/BBarBrickUP.h>
void *OPS_BBarBrickUP()
{
  if (OPS_GetNDM() != 3 || OPS_GetNDF() != 4) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element brickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? "
              "matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return 0;
  }

  // BBarBrickUPId, Nod[8], matID
  int tags[10];
  int num = 10;
  if (OPS_GetIntInput(&num, tags) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(tags[9]);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << tags[9];
    opserr << "\nBBarBrickUP element: " << tags[0] << endln;
    return 0;
  }

  // bk, r, perm1, perm2, perm3
  double data[5];
  num = 5;
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: invalid double input\n";
    return 0;
  }

  // b1, b2, b3
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

  return new BBarBrickUP(tags[0], tags[1], tags[2], tags[3], tags[4], tags[5],
                         tags[6], tags[7], tags[8], *mat, data[0], data[1],
                         data[2], data[3], data[4], opt[0], opt[1], opt[2]);
}
