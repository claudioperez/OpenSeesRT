

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/QzLiq1.h>
void *OPS_QzLiq1(G3_Runtime *rt)
{
  UniaxialMaterial *theMat = 0;

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial QzLiq1 tag? qzType? qult? z50? suction? "
              "dashpot? alpha? solidElem1? solidElem2?\n";
    opserr << "or: uniaxialMaterial QzLiq1 tag? qzType? qult? z50? suction? "
              "dashpot? alpha? -timeSeries seriesTag?\n";
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
    theMat = new QzLiq1(idata[0], idata[1], ddata[0], ddata[1], ddata[2],
                        ddata[3], ddata[4], theDomain, theSeries);

  } else {

    // back one arg
    OPS_ResetCurrentInputArg(-1);

    int eleTags[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, eleTags) < 0) {
      opserr << "WARNING invalid element tags\n";
      return 0;
    }

    theMat = new QzLiq1(idata[0], idata[1], ddata[0], ddata[1], ddata[2],
                        ddata[3], ddata[4], eleTags[0], eleTags[1], theDomain);
  }

  return theMat;
}
