

#include <g3_api.h>


#include <SRC/element/PFEMElement/PFEMElement2DFIC.h>
void *OPS_PFEMElement2DFIC()
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

  // rho, mu, b1, b2
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 6)
    numdata = 6;
  double data[6] = {0, 0, 0, 0, 1.0, 1e9};
  if (OPS_GetDoubleInput(&numdata, data) < 0)
    return 0;

  return new PFEMElement2DFIC(idata[0], idata[1], idata[2], idata[3], data[0],
                              data[1], data[2], data[3], data[4]);
}
