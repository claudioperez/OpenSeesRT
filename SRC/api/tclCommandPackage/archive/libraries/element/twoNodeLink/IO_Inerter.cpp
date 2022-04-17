

#include <g3_api.h>


#include <SRC/element/twoNodeLink/Inerter.h>
void *OPS_Inerter()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: inerter eleTag iNode jNode -dir dirs -inertance ib "
              "<-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-doRayleigh> "
              "<-damp cb> <-mass m>\n";
    return 0;
  }

  // tags
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  // dirs
  const char *type = OPS_GetString();
  if (strcmp(type, "-dir") != 0 && strcmp(type, "-dof") != 0) {
    opserr << "WARNING expecting -dir dirs\n";
    return 0;
  }
  ID dirs(ndf);
  int numDIR = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    int dir;
    numdata = 1;
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (OPS_GetIntInput(&numdata, &dir) < 0) {
      if (numArgs > OPS_GetNumRemainingInputArgs()) {
        // move current arg back by one
        OPS_ResetCurrentInputArg(-1);
      }
      break;
    }
    if (dir < 1 || ndf < dir) {
      opserr << "WARNING invalid direction ID\n";
      return 0;
    }
    dirs(numDIR++) = dir - 1;
  }
  dirs.resize(numDIR);

  // inertance matrix terms
  type = OPS_GetString();
  if (strcmp(type, "-inertance") != 0 && strcmp(type, "-inertia") != 0) {
    opserr << "WARNING expecting -inertance ib\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < numDIR * numDIR) {
    opserr << "WARNING wrong number of ib values specified\n";
    return 0;
  }
  numdata = 1;
  Matrix ib(numDIR, numDIR);
  for (int i = 0; i < numDIR; i++) {
    for (int j = 0; j < numDIR; j++) {
      if (OPS_GetDoubleInput(&numdata, &ib(i, j)) < 0) {
        opserr << "WARNING invalid inertance value\n";
        return 0;
      }
    }
  }

  // options
  Vector x, y, Mratio;
  int doRayleigh = 0;
  Matrix *cb = 0;
  double mass = 0.0;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    return new Inerter(idata[0], ndm, idata[1], idata[2], dirs, ib);
  }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    type = OPS_GetString();
    if (strcmp(type, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING: insufficient arguments after -orient\n";
        return 0;
      }
      numdata = 3;
      x.resize(3);
      if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
        opserr << "WARNING: invalid -orient values\n";
        return 0;
      }
      if (OPS_GetNumRemainingInputArgs() < 3) {
        y = x;
        x = Vector();
        continue;
      }
      y.resize(3);
      if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
        y = x;
        x = Vector();
        continue;
      }
    } else if (strcmp(type, "-pDelta") == 0) {
      Mratio.resize(4);
      Mratio.Zero();
      numdata = 4;
      double *ptr = &Mratio(0);
      if (ndm == 2) {
        numdata = 2;
        ptr += 2;
      }
      if (OPS_GetNumRemainingInputArgs() < numdata) {
        opserr << "WARNING: insufficient data for -pDelta\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&numdata, ptr) < 0) {
        opserr << "WARNING: invalid -pDelta value\n";
        return 0;
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(type, "-damp") == 0) {
      if (OPS_GetNumRemainingInputArgs() < numDIR * numDIR) {
        opserr << "WARNING wrong number of cb values specified\n";
        return 0;
      }
      numdata = 1;
      double cij;
      cb = new Matrix(numDIR, numDIR);
      for (int i = 0; i < numDIR; i++) {
        for (int j = 0; j < numDIR; j++) {
          if (OPS_GetDoubleInput(&numdata, &cij) < 0) {
            opserr << "WARNING invalid damping value\n";
            delete cb;
            return 0;
          }
          (*cb)(i, j) = cij;
        }
      }
    } else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WANRING: insufficient mass value\n";
        return 0;
      }
      numdata = 1;
      if (OPS_GetDoubleInput(&numdata, &mass) < 0) {
        opserr << "WANRING: invalid -mass value\n";
        return 0;
      }
    }
  }

  // create object
  Element *theEle = new Inerter(idata[0], ndm, idata[1], idata[2], dirs, ib, y,
                                x, Mratio, doRayleigh, cb, mass);

  // clean up memory
  delete cb;

  return theEle;
}
