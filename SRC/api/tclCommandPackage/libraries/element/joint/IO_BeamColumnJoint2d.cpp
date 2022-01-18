

#include <g3_api.h>


#include <SRC/element/joint/BeamColumnJoint2d.h>
void *OPS_BeamColumnJoint2d()
{
  if (OPS_GetNumRemainingInputArgs() < 18) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element beamColumnJoint eleTag? node1? node2? node3? "
              "node4? matTag1? matTag2? matTag3?\n";
    opserr << "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? matTag10? "
              "matTag11? matTag12? matTag13?\n";
    opserr << "<ElementHeightFactor? ElementWidthFactor?>\n";
    return 0;
  }

  int idata[18];
  int numdata = 18;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  double data[2] = {1.0, 1.0};
  numdata = 2;
  if (OPS_GetNumRemainingInputArgs() > 1) {
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
      opserr << "WARNING: invalid double inputs\n";
      return 0;
    }
  }

  UniaxialMaterial *mats[13];
  for (int i = 0; i < 13; i++) {
    mats[i] = OPS_getUniaxialMaterial(idata[5 + i]);
    if (mats[i] == 0) {
      opserr << "WARNING: material " << idata[5 + i] << " is not defined\n";
      return 0;
    }
  }

  return new BeamColumnJoint2d(
      idata[0], idata[1], idata[2], idata[3], idata[4], *mats[0], *mats[1],
      *mats[2], *mats[3], *mats[4], *mats[5], *mats[6], *mats[7], *mats[8],
      *mats[9], *mats[10], *mats[11], *mats[12], data[0], data[1]);
}
