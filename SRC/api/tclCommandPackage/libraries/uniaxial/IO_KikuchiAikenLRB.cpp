

#include <g3_api.h>


#include <SRC/material/uniaxial/KikuchiAikenLRB.h>
void *OPS_KikuchiAikenLRB(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 9) {
    opserr << "WARNING invalid number of arguments\n";
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid KikuchiAikenHDR tag\n";
    return 0;
  }

  double ddata[7];
  numdata = 7;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double temp = 15.0;
  double ddata2[2] = {1, 1};
  double ddata3[2] = {1, 1};

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-coKQ") == 0) {
      if (OPS_GetNumRemainingInputArgs() >= 2) {
        numdata = 2;
        if (OPS_GetDoubleInput(&numdata, ddata2) < 0) {
          opserr << "WARNING invalid double inputs\n";
          return 0;
        }
      }
    } else if (strcmp(opt, "-coMSS") == 0) {
      if (OPS_GetNumRemainingInputArgs() >= 2) {
        numdata = 2;
        if (OPS_GetDoubleInput(&numdata, ddata3) < 0) {
          opserr << "WARNING invalid double inputs\n";
          return 0;
        }
      }
    } else if (strcmp(opt, "-T") == 0) {
      if (OPS_GetNumRemainingInputArgs() >= 1) {
        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &temp) < 0) {
          opserr << "WARNING invalid temp\n";
          return 0;
        }
      }

    } else {
      opserr << "WARNING invalid optional arguments \n";
      return 0;
    }
  }

  for (int i = 0; i < 2; i++) {
    if (ddata2[i] == 0.0)
      ddata2[i] = 1.0;
  }
  for (int i = 0; i < 2; i++) {
    if (ddata3[i] == 0.0)
      ddata3[i] = 1.0;
  }

  return new KikuchiAikenLRB(idata[0], idata[1], ddata[0], ddata[1], ddata[2],
                             ddata[3], ddata[4], ddata[5], ddata[6], temp,
                             ddata2[0], ddata2[1], ddata3[0], ddata3[1]);
}
