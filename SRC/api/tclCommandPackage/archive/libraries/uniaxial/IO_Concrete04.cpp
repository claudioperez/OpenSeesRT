

#include <g3_api.h>


#include <SRC/material/uniaxial/Concrete04.h>
void *OPS_Concrete04(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? Ec0?";
    opserr << " <ft? etu? <beta?> >\n";
    return 0;
  }

  int type = 1;

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double data\n";
    return 0;
  }

  numdata = OPS_GetNumRemainingInputArgs();
  double data2[2];
  if (numdata > 1) {
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, data2) < 0) {
      opserr << "WARNING invalid double data\n";
      return 0;
    }
    type = 2;
  }

  numdata = OPS_GetNumRemainingInputArgs();
  double beta;
  if (numdata > 0) {
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &beta)) {
      opserr << "WARNING invalid double data\n";
      return 0;
    }
    type = 3;
  }

  UniaxialMaterial *mat = 0;
  if (type == 1) {
    mat = new Concrete04(tag, data[0], data[1], data[2], data[3]);
  } else if (type == 2) {
    mat = new Concrete04(tag, data[0], data[1], data[2], data[3], data2[0],
                         data2[1]);
  } else if (type == 3) {
    mat = new Concrete04(tag, data[0], data[1], data[2], data[3], data2[0],
                         data2[1], beta);
  }

  if (mat == 0) {
    opserr << "WARNING: failed to create Concrete04 material\n";
    return 0;
  }

  return mat;
}
