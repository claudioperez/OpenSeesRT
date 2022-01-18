

#include <g3_api.h>


#include <SRC/element/PFEMElement/TaylorHood2D.h>
void *OPS_TaylorHood2D()
{
  int num = OPS_GetNumRemainingInputArgs();
  if (num < 4) {
    opserr << "WARNING: insufficient number of arguments -- TaylorHood2D:";
    opserr << "tag, nd1, nd2, nd3, (rho, mu, b1, b2, thk, kappa)\n";
    return 0;
  }

  // tag, nd1, nd2, nd3
  num = 4;
  int idata[4];
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: failed to read integers -- TaylorHood2D\n";
    return 0;
  }

  // (rho, mu, b1, b2, thinkness, kappa)
  num = OPS_GetNumRemainingInputArgs();
  if (num > 6)
    num = 6;
  double data[6] = {1000.0, 1e-3, 0, -9.81, 1.0, 2.15e9};
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: failed to read doubles -- TaylorHood2D\n";
    return 0;
  }

  return new TaylorHood2D(idata[0], idata[1], idata[2], idata[3], data[0],
                          data[1], data[2], data[3], data[4], data[5]);
}
