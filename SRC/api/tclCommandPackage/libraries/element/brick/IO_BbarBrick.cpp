

#include <g3_api.h>


#include <SRC/element/brick/BbarBrick.h>
void *OPS_BbarBrick()
{
  if (OPS_GetNumRemainingInputArgs() < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element Brick eleTag? Node1? Node2? Node3? Node4? Node5? "
              "Node6? Node7? Node 8? matTag?\n";
    return 0;
  }

  int idata[10];
  int num = 10;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[9]);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << idata[9];
    opserr << "\nBrick element: " << idata[0] << endln;
  }

  double data[3] = {0, 0, 0};
  num = OPS_GetNumRemainingInputArgs();
  if (num > 3) {
    num = 3;
  }
  if (num > 0) {
    if (OPS_GetDoubleInput(&num, data) < 0) {
      opserr << "WARNING: invalid double data\n";
      return 0;
    }
  }

  return new BbarBrick(idata[0], idata[1], idata[2], idata[3], idata[4],
                       idata[5], idata[6], idata[7], idata[8], *mat, data[0],
                       data[1], data[2]);
}
