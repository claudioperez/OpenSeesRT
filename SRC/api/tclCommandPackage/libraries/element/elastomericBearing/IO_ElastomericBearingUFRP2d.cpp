

#include <g3_api.h>


#include <SRC/element/elastomericBearing/ElastomericBearingUFRP2d.h>
void *OPS_ElastomericBearingUFRP2d()
{
  int ndf = OPS_GetNDF();

  // check plane frame problem has 3 dof per node
  if (ndf != 3) {
    opserr << "WARNING invalid ndf: " << ndf;
    opserr << ", for plane problem need 3 - elastomericBearingUFRP\n";
    return 0;
  }

  // check the number of arguments is correct
  if (OPS_GetNumRemainingInputArgs() < 18) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: elastomericBearingUFRP eleTag iNode jNode uy a1 a2 a3 a4 "
              "a5 b c eta beta gamma -P matTag -Mz matTag <-orient x1 x2 x3 y1 "
              "y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter "
              "maxIter tol>\n";
    return 0;
  }

  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs \n";
    return 0;
  }
  int tag = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];

  double data[11];
  numdata = 1;
  ;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs \n";
    return 0;
  }
  double uy = data[0];
  double a1 = data[1];
  double a2 = data[2];
  double a3 = data[3];
  double a4 = data[4];
  double a5 = data[5];
  double b = data[6];
  double c = data[7];
  double eta = data[8];
  double beta = data[9];
  double gamma = data[10];

  UniaxialMaterial *theMaterials[2] = {0, 0};
  const char *flag = OPS_GetString();
  if (strcmp(flag, "-P") != 0) {
    opserr << "WARNING -P is expected\n";
    return 0;
  }
  int matTag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid matTag\n";
    return 0;
  }
  theMaterials[0] = OPS_getUniaxialMaterial(matTag);
  if (theMaterials[0] == 0) {
    opserr << "WARNING material model not found\n";
    opserr << "uniaxialMaterial: " << matTag << endln;
    opserr << "elastomericBearingUFRP element: " << tag << endln;
    return 0;
  }

  flag = OPS_GetString();
  if (strcmp(flag, "-Mz") != 0) {
    opserr << "WARNING -Mz is expected\n";
    return 0;
  }
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid matTag\n";
    return 0;
  }
  theMaterials[1] = OPS_getUniaxialMaterial(matTag);
  if (theMaterials[1] == 0) {
    opserr << "WARNING material model not found\n";
    opserr << "uniaxialMaterial: " << matTag << endln;
    opserr << "elastomericBearingUFRP element: " << tag << endln;
    return 0;
  }

  // get the id and end nodes
  double shearDistI = 0.5;
  int doRayleigh = 0;
  double mass = 0.0;
  int maxIter = 25;
  double tol = 1E-12;

  // check for optional arguments
  Vector x;
  Vector y;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    flag = OPS_GetString();
    if (strcmp(flag, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING insufficient arguments after -orient flag\n";
        opserr << "elastomericBearingUFRP element: " << tag << endln;
        return 0;
      }
      x.resize(3);
      y.resize(3);
      numdata = 3;
      if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
        opserr << "WARNING invalid -orient value\n";
        opserr << "elastomericBearingUFRP element: " << tag << endln;
        return 0;
      }
      if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
        opserr << "WARNING invalid -orient value\n";
        opserr << "elastomericBearingUFRP element: " << tag << endln;
        return 0;
      }
    } else if (strcmp(flag, "-shearDist") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numdata, &shearDistI) < 0) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "elastomericBearingUFRP element: " << tag << endln;
          return 0;
        }
      }
    } else if (strcmp(flag, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(flag, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numdata, &mass) < 0) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "elastomericBearingUFRP element: " << tag << endln;
          return 0;
        }
      }
    } else if (strcmp(flag, "-massiter") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 1) {
        if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
          opserr << "WARNING invalid maxIter value\n";
          opserr << "elastomericBearingUFRP element: " << tag << endln;
          return 0;
        }
        if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
          opserr << "WARNING invalid tol value\n";
          opserr << "elastomericBearingUFRP element: " << tag << endln;
          return 0;
        }
      }
    }
  }

  // now create the elastomericBearingUFRP
  return new ElastomericBearingUFRP2d(
      tag, iNode, jNode, uy, a1, a2, a3, a4, a5, b, c, theMaterials, y, x, eta,
      beta, gamma, shearDistI, doRayleigh, mass, maxIter, tol);
}
