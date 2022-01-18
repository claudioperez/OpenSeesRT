

#include <g3_api.h>


#include <SRC/element/elastomericBearing/ElastomericBearingBoucWen2d.h>
void *OPS_ElastomericBearingBoucWen2d()
{
  int ndf = OPS_GetNDF();
  if (ndf != 3) {
    opserr << "WARNING invalid ndf: " << ndf;
    opserr << ", for plane problem need 3 - elastomericBearingBoucWen\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: elastomericBearingBoucWen eleTag iNode jNode kInit qd "
              "alpha1 alpha2 mu eta beta gamma -P matTag -Mz matTag <-orient "
              "x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> "
              "<-iter maxIter tol>\n";
    return 0;
  }

  // tags
  int idata[3];
  int num = 3;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  // data
  double data[8];
  num = 8;
  if (OPS_GetDoubleInput(&num, data) < 0) {
    opserr << "WARNING: invalid double inputs\n";
    return 0;
  }

  // materials
  UniaxialMaterial *mats[2] = {0, 0};
  const char *type = OPS_GetString();
  if (strcmp(type, "-P") != 0) {
    opserr << "WARNING: want -P\n";
    return 0;
  }
  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING: invalid matTag\n";
    return 0;
  }
  mats[0] = OPS_getUniaxialMaterial(matTag);
  if (mats[0] == 0) {
    opserr << "WARNING: material not found\n";
    return 0;
  }

  type = OPS_GetString();
  if (strcmp(type, "-Mz") != 0) {
    opserr << "WARNING: want -Mz\n";
    return 0;
  }
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING: invalid matTag\n";
    return 0;
  }
  mats[1] = OPS_getUniaxialMaterial(matTag);
  if (mats[1] == 0) {
    opserr << "WARNING: material not found\n";
    return 0;
  }

  // options
  Vector x, y;
  double sDistI = 0.5;
  int doRayleigh = 0;
  double mass = 0.0;
  int maxIter = 25;
  double tol = 1e-12;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    return new ElastomericBearingBoucWen2d(
        idata[0], idata[1], idata[2], data[0], data[1], data[2], mats, y, x,
        data[3], data[4], data[5], data[6], data[7], sDistI, doRayleigh, mass,
        maxIter, tol);
  }
  while (OPS_GetNumRemainingInputArgs() > 0) {
    type = OPS_GetString();
    if (strcmp(type, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING: insufficient arguments after -orient\n";
        return 0;
      }
      num = 3;
      x.resize(3);
      if (OPS_GetDoubleInput(&num, &x(0)) < 0) {
        opserr << "WARNING: invalid orient value\n";
        return 0;
      }
      y.resize(3);
      if (OPS_GetDoubleInput(&num, &y(0)) < 0) {
        opserr << "WARNING: invalid orient value\n";
        return 0;
      }
    } else if (strcmp(type, "-shearDist") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: insufficient args\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &sDistI) < 0) {
        opserr << "WARNING: invalid shearDist\n";
        return 0;
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: insufficient args\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &mass) < 0) {
        opserr << "WARNING: invalid mass\n";
        return 0;
      }
    } else if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING: insufficient args\n";
        return 0;
      }
      num = 1;
      if (OPS_GetIntInput(&num, &maxIter) < 0) {
        opserr << "WARNING: invalid maxIter\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&num, &tol) < 0) {
        opserr << "WARNING: invalid tol\n";
        return 0;
      }
    }
  }

  return new ElastomericBearingBoucWen2d(
      idata[0], idata[1], idata[2], data[0], data[1], data[2], mats, y, x,
      data[3], data[4], data[5], data[6], data[7], sDistI, doRayleigh, mass,
      maxIter, tol);
}
