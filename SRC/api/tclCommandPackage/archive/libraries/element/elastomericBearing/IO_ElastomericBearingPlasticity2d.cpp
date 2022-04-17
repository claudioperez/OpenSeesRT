

#include <g3_api.h>


#include <SRC/element/elastomericBearing/ElastomericBearingPlasticity2d.h>
void *OPS_ElastomericBearingPlasticity2d()
{
  int ndf = OPS_GetNDF();
  if (ndf != 3) {
    opserr << "WARNING invalid ndf: " << ndf;
    opserr << ", for plane problem need 3 - elastomericBearing\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 12) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: elastomericBearing eleTag iNode jNode kInit qd alpha1 "
              "alpha2 mu -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> "
              "<-shearDist sDratio> <-doRayleigh> <-mass m>\n";
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
  double data[5];
  num = 5;
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
  Vector x, y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;
  double sDistI = 0.5;
  int doRayleigh = 0;
  double mass = 0.0;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    return new ElastomericBearingPlasticity2d(
        idata[0], idata[1], idata[2], data[0], data[1], data[2], mats, y, x,
        data[3], data[4], sDistI, doRayleigh, mass);
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
    }
  }

  return new ElastomericBearingPlasticity2d(
      idata[0], idata[1], idata[2], data[0], data[1], data[2], mats, y, x,
      data[3], data[4], sDistI, doRayleigh, mass);
}
