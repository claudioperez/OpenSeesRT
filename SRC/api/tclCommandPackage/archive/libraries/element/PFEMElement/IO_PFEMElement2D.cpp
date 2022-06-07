

#include <g3_api.h>


#include <SRC/element/PFEMElement/PFEMElement2D.h>
void *OPS_PFEMElement2D()
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 8) {
    opserr << "WARNING: insufficient number of arguments\n";
    return 0;
  }

  // tag, nd1, nd2, nd3
  numdata = 4;
  int idata[4];
  if (OPS_GetIntInput(&numdata, idata) < 0)
    return 0;

  // rho, mu, b1, b2, (thickness, kappa, lumped)
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 7)
    numdata = 7;
  double data[7] = {0, 0, 0, 0, 1.0, -1, 1};
  if (OPS_GetDoubleInput(&numdata, data) < 0)
    return 0;

  return new PFEMElement2D(idata[0], idata[1], idata[2], idata[3], data[0],
                           data[1], data[2], data[3], data[4], data[5],
                           data[6]);
}
