

#include <g3_api.h>


#include <SRC/material/uniaxial/AxialSpHD.h>
void *OPS_AxialSpHD(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING invalid number of arguments\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid AxialSp tag\n";
    return 0;
  }

  double data[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double opt[6] = {1, 1, 1, 1, 0, 1};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 6)
    numdata = 6;
  if (OPS_GetDoubleInput(&numdata, opt) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  return new AxialSpHD(tag, data[0], data[1], data[2], opt[0], opt[1], opt[2],
                       opt[3], opt[4], opt[5]);
}
