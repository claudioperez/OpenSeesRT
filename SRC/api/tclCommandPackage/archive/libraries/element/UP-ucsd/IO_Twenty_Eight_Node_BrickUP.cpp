

#include <g3_api.h>


#include <SRC/element/UP-ucsd/Twenty_Eight_Node_BrickUP.h>
void *OPS_TwentyEightNodeBrickUP()
{
  if (OPS_GetNDM() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with 20_8_BrickUP element\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 27) {
    opserr << "WARNING insufficient arguments\n";
    opserr
        << "Want: element 20_8_BrickUP eleTag? Node1? ... Node20? thk? type? "
           "matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return 0;
  }

  // brickUPId, Node[20], matID
  int tags[22];
  int num = 22;
  if (OPS_GetIntInput(&num, tags) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(tags[21]);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << tags[21];
    opserr << "\nBrick element: " << tags[0] << endln;
  }

  // bk, r, perm1, perm2, perm3
  double data[5];
  num = 5;
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: invalid double input\n";
    return 0;
  }

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

  return new TwentyEightNodeBrickUP(
      tags[0], tags[1], tags[2], tags[3], tags[4], tags[5], tags[6], tags[7],
      tags[8], tags[9], tags[10], tags[11], tags[12], tags[13], tags[14],
      tags[15], tags[16], tags[17], tags[18], tags[19], tags[20], *mat, data[0],
      data[1], data[2], data[3], data[4], opt[0], opt[1], opt[2]);
}
