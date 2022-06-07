

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/TzLiq1.h>
void *OPS_TzLiq1(G3_Runtime *rt)
{
  UniaxialMaterial *theMat = 0;

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial TzLiq1 tag? tzType? tult? z50? dashpot? "
              "solidElem1? solidElem2?\n";
    opserr << "or: uniaxialMaterial TzLiq1 tag? tzType? tult? z50? dashpot? "
              "-timeSeries seriesTag?\n";
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  double ddata[3];
  numdata = 3;
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
    theMat = new TzLiq1(idata[0], MAT_TAG_TzLiq1, idata[1], ddata[0], ddata[1],
                        ddata[2], theDomain, theSeries);

  } else {

    // back one arg
    OPS_ResetCurrentInputArg(-1);

    int eleTags[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, eleTags) < 0) {
      opserr << "WARNING invalid element tags\n";
      return 0;
    }

    theMat = new TzLiq1(idata[0], MAT_TAG_TzLiq1, idata[1], ddata[0], ddata[1],
                        ddata[2], eleTags[0], eleTags[1], theDomain);
  }

  return theMat;
}
