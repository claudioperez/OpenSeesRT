

#include <g3_api.h>


#include <SRC/element/PFEMElement/MINI.h>
void *OPS_MINI()
{
  int ndm = OPS_GetNDM();
  int num = OPS_GetNumRemainingInputArgs();
  if (ndm == 2) {
    if (num < 10) {
      opserr << "WARNING: insufficient number of arguments -- MINI:";
      opserr << "tag, nd1, nd2, nd3, rho, mu, b1, b2, thk, kappa\n";
      return 0;
    }

  } else if (ndm == 3) {
    if (num < 12) {
      opserr << "WARNING: insufficient number of arguments -- MINI:";
      opserr << "tag, nd1, nd2, nd3, nd4, rho, mu, b1, b2, b3, thk, kappa\n";
      return 0;
    }
  }

  // tag, nd1, nd2, nd3
  if (ndm == 2) {
    num = 4;
  } else if (ndm == 3) {
    num = 5;
  }
  int idata[5];
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: failed to read integers -- MINI\n";
    return 0;
  }

  // (rho, mu, b1, b2,b3 thinkness, kappa)
  if (ndm == 2) {
    num = 6;
  } else if (ndm == 3) {
    num = 7;
  }
  double data[7];
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: failed to read doubles -- MINI\n";
    return 0;
  }

  if (ndm == 2) {

    return new MINI(idata[0], idata[1], idata[2], idata[3], data[0], data[1],
                    data[2], data[3], data[4], data[5]);
  } else if (ndm == 3) {
    return new MINI(idata[0], idata[1], idata[2], idata[3], idata[4], data[0],
                    data[1], data[2], data[3], data[4], data[5], data[6]);
  }

  return 0;
}
