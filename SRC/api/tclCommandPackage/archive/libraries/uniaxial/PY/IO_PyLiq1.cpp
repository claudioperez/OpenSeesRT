

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/PyLiq1.h>
void *OPS_PyLiq1(G3_Runtime *rt)
{
  UniaxialMaterial *theMat = 0;

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? "
              "dashpot? pRes? solidElem1? solidElem2?\n";
    opserr << "or: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? "
              "dashpot? -timeSeries seriesTag?\n";
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  double ddata[5];
  numdata = 5;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  const char *arg = OPS_GetString();
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;

  if (strcmp(arg, "-timeSeries") == 0) {
    int tsTag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
      opserr << "WARNING invalid time series tag\n";
      return 0;
    }
    TimeSeries *theSeries = OPS_getTimeSeries(tsTag);
    theMat = new PyLiq1(idata[0], MAT_TAG_PyLiq1, idata[1], ddata[0], ddata[1],
                        ddata[2], ddata[3], ddata[4], theDomain, theSeries);

  } else {

    // back one arg
    OPS_ResetCurrentInputArg(-1);

    int eleTags[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, eleTags) < 0) {
      opserr << "WARNING invalid element tags\n";
      return 0;
    }

    theMat = new PyLiq1(idata[0], MAT_TAG_PyLiq1, idata[1], ddata[0], ddata[1],
                        ddata[2], ddata[3], ddata[4], eleTags[0], eleTags[1],
                        theDomain);
  }

  return theMat;
}
