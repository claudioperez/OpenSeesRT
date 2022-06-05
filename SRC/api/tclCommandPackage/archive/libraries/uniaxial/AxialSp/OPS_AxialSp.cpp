

#include <g3_api.h>


#include <SRC/material/uniaxial/AxialSp.h>
void *OPS_AxialSp(G3_Runtime *rt)
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

  double opt[4] = {0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 4)
    numdata = 4;
  if (OPS_GetDoubleInput(&numdata, opt) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  return new AxialSp(tag, data[0], data[1], data[2], opt[0], opt[1], opt[2],
                     opt[3]);
}
