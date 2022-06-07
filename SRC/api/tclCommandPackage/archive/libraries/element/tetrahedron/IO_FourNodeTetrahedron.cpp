

#include <g3_api.h>


#include <SRC/element/tetrahedron/FourNodeTetrahedron.h>
void *OPS_FourNodeTetrahedron()
{
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element FourNodeTetrahedron eleTag? Node1? Node2? Node3? "
              "Node4? matTag?\n";
    return 0;
  }

  int idata[6];
  int num = 6;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[5]);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << idata[9];
    opserr << "\nFourNodeTetrahedron element: " << idata[0] << endln;
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
  return new FourNodeTetrahedron(idata[0], idata[1], idata[2], idata[3],
                                 idata[4], *mat, data[0], data[1], data[2]);
}
