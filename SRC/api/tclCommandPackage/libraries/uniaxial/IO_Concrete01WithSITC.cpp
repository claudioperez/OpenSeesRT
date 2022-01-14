

#include <g3_api.h>


#include <SRC/material/uniaxial/Concrete01WithSITC.h>
void *OPS_Concrete01WithSITC(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Concrete01WithSITC tag? ";
    opserr << "fpc? epsc0? fpcu? epscu? <endStrainSITC?>\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, data)) {
    opserr << "WARNING invalid double data\n";
    return 0;
  }
  UniaxialMaterial *mat = 0;

  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 0) {
    double endStrainSITC;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &endStrainSITC) < 0) {
      opserr << "WARNING invalid double data\n";
      return 0;
    }
    mat = new Concrete01WithSITC(tag, data[0], data[1], data[2], data[3],
                                 endStrainSITC);
  } else {
    mat = new Concrete01WithSITC(tag, data[0], data[1], data[2], data[3]);
  }

  if (mat == 0) {
    opserr << "WARNING: failed to create Concrete01WithSITC material\n";
    return 0;
  }

  return mat;
}
