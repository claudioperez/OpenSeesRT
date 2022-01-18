

#include <g3_api.h>


#include <SRC/element/brick/Twenty_Node_Brick.h>
void *OPS_Twenty_Node_Brick()
{
  if (OPS_GetNDM() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with 20NodeBrick element\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 22) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element 20NodeBrick eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? "
              "matTag? <b1? b2? b3?>\n";
    return 0;
  }

  // brickId, Nod[20], matID
  int idata[22];
  int num = 22;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[21]);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << idata[21];
    opserr << "\nBrick element: " << idata[0] << endln;
  }

  // b1, b2, b3
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

  return new Twenty_Node_Brick(
      idata[0], idata[1], idata[2], idata[3], idata[4], idata[5], idata[6],
      idata[7], idata[8], idata[9], idata[10], idata[11], idata[12], idata[13],
      idata[14], idata[15], idata[16], idata[17], idata[18], idata[19],
      idata[20], *mat, data[0], data[1], data[2]);
}
