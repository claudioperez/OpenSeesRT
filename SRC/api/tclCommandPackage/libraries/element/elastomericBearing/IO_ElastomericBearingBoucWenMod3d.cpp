

#include <g3_api.h>


#include <SRC/element/elastomericBearing/ElastomericBearingBoucWenMod3d.h>
void *OPS_ElastomericBearingBoucWenMod3d(void)
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 18) {
    opserr
        << "Want: ElastomericBearingBoucWenMod eleTag iNode jNode kInit fy Gr "
           "Kbulk D1 D2 ts tr n alpha1 alpha2 mu eta beta gamma <-PMod a1 a2 > "
           "<-TMod T b1 b2 b3 b4> <-shearDist sDratio> <-doRayleigh> <-mass m> "
           "<-iter maxIter tol> <-orient <x1 x2 x3> y1 y2 y3>\n";
    return 0;
  }
  int tag;
  int iNode, jNode;
  double kInit, fy, alpha1, Gr, Kbulk, D1, D2, ts, tr;
  int n;

  double alpha2 = 0.0;
  double mu = 2.0;
  double eta = 1.0;
  double beta = 0.5;
  double gamma = 0.5;
  double a1 = 0.0;
  double a2 = 1.0;
  double T = 23.0;
  double b1 = 1.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double shearDistI = 0.5;

  int doRayleigh = 0;
  double mass = 0.0;
  int maxIter = 25;
  double tol = 1E-12;

  Vector x(0);
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;

  int iData[3];
  double dData[15];

  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d error reading integer "
              "data\n";
    return 0;
  }
  tag = iData[0];
  iNode = iData[1];
  jNode = iData[2];

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
           << " error reading data\n";
    return 0;
  }
  kInit = dData[0];
  fy = dData[1];
  Gr = dData[2];
  Kbulk = dData[3];
  D1 = dData[4];
  D2 = dData[5];
  ts = dData[6];
  tr = dData[7];

  numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
           << " error reading data\n";
    return 0;
  }
  n = iData[0];

  numData = 6;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
           << " error reading data\n";
    return 0;
  }
  alpha1 = dData[0];
  alpha2 = dData[1];
  mu = dData[2];
  eta = dData[3];
  beta = dData[4];
  gamma = dData[5];

  // while not done, loop over options
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *option = OPS_GetString();
    if (strcmp(option, "-PMod") == 0) {
      numData = 2;
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -PMod data\n";
        break;
      }
      a1 = dData[0];
      a2 = dData[1];
    } else if (strcmp(option, "-TMod") == 0) {
      numData = 4;
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -TMod data\n";
        break;
      }
      b1 = dData[0];
      b1 = dData[1];
      b1 = dData[2];
      b1 = dData[3];
    } else if (strcmp(option, "-shearDist") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -shearDist data\n";
        break;
      }
      shearDistI = dData[0];
    } else if (strcmp(option, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(option, "-mass") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -mass data\n";
        break;
      }
      mass = dData[0];
    } else if (strcmp(option, "-iter") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -iter data\n";
        break;
      }
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -iter data\n";
        break;
      }
      maxIter = iData[0];
      tol = dData[0];
    } else if (strcmp(option, "-orient") == 0) {
      int numOrient = OPS_GetNumRemainingInputArgs();
      numData = numOrient;
      if (numData != 3 || numData != 6) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -orient data, need 3 or 6 values\n";
      }
      if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag
               << " error reading -orient data\n";
        break;
      }
      if (numOrient == 3) {
        y(0) = dData[0];
        y(1) = dData[1];
        y(2) = dData[2];
      } else if (numOrient == 6) {
        x.resize(3);
        x(0) = dData[0];
        x(1) = dData[1];
        x(2) = dData[2];
        y(0) = dData[3];
        y(1) = dData[4];
        y(2) = dData[5];
      }
    }
  }
  theElement = new ElastomericBearingBoucWenMod3d(
      tag, iNode, jNode, kInit, fy, Gr, Kbulk, D1, D2, ts, tr, n, alpha1,
      alpha2, mu, eta, beta, gamma, a1, a2, T, b1, b2, b3, b4, y, x, shearDistI,
      doRayleigh, mass, maxIter, tol);

  return theElement;
}
