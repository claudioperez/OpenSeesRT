

#include <g3_api.h>


#include <SRC/element/twoNodeLink/TwoNodeLink.h>
void *OPS_TwoNodeLink()
{
  int ndm = OPS_GetNDM();
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: twoNodeLink eleTag iNode jNode -mat matTags -dir dirs "
              "<-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-shearDist "
              "sDratios> <-doRayleigh> <-mass m>\n";
    return 0;
  }

  // tags
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  // mats
  const char *type = OPS_GetString();
  if (strcmp(type, "-mat") != 0) {
    opserr << "WARNING expecting -mat matTags\n";
    return 0;
  }
  std::vector<UniaxialMaterial *> mats;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    int mattag;
    numdata = 1;
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (OPS_GetIntInput(&numdata, &mattag) < 0) {
      if (numArgs > OPS_GetNumRemainingInputArgs()) {
        // move current arg back by one
        OPS_ResetCurrentInputArg(-1);
      }
      break;
    }
    UniaxialMaterial *mat = OPS_getUniaxialMaterial(mattag);
    if (mat == 0) {
      opserr << "WARNING material model not found\n";
      opserr << "uniaxialMaterial " << mattag << endln;
      return 0;
    }
    mats.push_back(mat);
  }

  // dirs
  type = OPS_GetString();
  if (strcmp(type, "-dir") != 0 && strcmp(type, "-dof") != 0) {
    opserr << "WARNING expecting -dir dirs\n";
    return 0;
  }
  ID dirs(int(mats.size()));
  if (OPS_GetNumRemainingInputArgs() < dirs.Size()) {
    opserr << "WARNING wrong number of directions specified\n";
    return 0;
  }
  numdata = dirs.Size();
  if (OPS_GetIntInput(&numdata, &dirs(0)) < 0) {
    opserr << "WARNING invalid direction ID\n";
    return 0;
  }
  for (int i = 0; i < numdata; i++)
    dirs(i)--;

  // options
  Vector x, y, Mratio, sDistI;
  int doRayleigh = 0;
  double mass = 0.0;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    return new TwoNodeLink(idata[0], ndm, idata[1], idata[2], dirs, &mats[0]);
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
    } else if (strcmp(type, "-shearDist") == 0) {
      sDistI.resize(2);
      numdata = 2;
      if (ndm == 2) {
        numdata = 1;
        sDistI(1) = 0.5;
      }
      if (OPS_GetNumRemainingInputArgs() < numdata) {
        opserr << "WARNING: insufficient data for -shearDist\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&numdata, &sDistI(0)) < 0) {
        opserr << "WARNING: invalid -shearDist value\n";
        return 0;
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      doRayleigh = 1;
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
  return new TwoNodeLink(idata[0], ndm, idata[1], idata[2], dirs, &mats[0], y,
                         x, Mratio, sDistI, doRayleigh, mass);
}
