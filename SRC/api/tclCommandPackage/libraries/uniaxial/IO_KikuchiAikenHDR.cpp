

#include <g3_api.h>


#include <SRC/material/uniaxial/KikuchiAikenHDR.h>
void *OPS_KikuchiAikenHDR(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING invalid number of arguments\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid KikuchiAikenHDR tag\n";
    return 0;
  }

  const char *arg = OPS_GetString();
  int tp;
  if ((strcmp(arg, "X0.6") == 0) || (strcmp(arg, "1") == 0)) {
    tp = 1;
  } else if ((strcmp(arg, "X0.6-0MPa") == 0) || (strcmp(arg, "2") == 0)) {
    tp = 2;
  } else if ((strcmp(arg, "X0.4") == 0) || (strcmp(arg, "3") == 0)) {
    tp = 3;
  } else if ((strcmp(arg, "X0.4-0MPa") == 0) || (strcmp(arg, "4") == 0)) {
    tp = 4;
  } else if ((strcmp(arg, "X0.3") == 0) || (strcmp(arg, "5") == 0)) {
    tp = 5;
  } else if ((strcmp(arg, "X0.3-0MPa") == 0) || (strcmp(arg, "6") == 0)) {
    tp = 6;
  } else {
    opserr << "WARNING invalid KikuchiAikenHDR tp\n";
    return 0;
  }

  double ddata[2];
  numdata = 2;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double ddata2[3] = {1, 1, 1};
  double ddata3[2] = {1, 1};

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-coGHU") == 0) {
      if (OPS_GetNumRemainingInputArgs() >= 3) {
        numdata = 3;
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
    } else {
      opserr << "WARNING invalid optional arguments \n";
      return 0;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (ddata2[i] == 0.0)
      ddata2[i] = 1.0;
  }
  for (int i = 0; i < 2; i++) {
    if (ddata3[i] == 0.0)
      ddata3[i] = 1.0;
  }

  return new KikuchiAikenHDR(tag, tp, ddata[0], ddata[1], ddata2[0], ddata2[1],
                             ddata2[2], ddata3[0], ddata3[1]);
}
