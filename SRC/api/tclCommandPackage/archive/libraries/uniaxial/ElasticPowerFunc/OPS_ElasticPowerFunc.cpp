

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticPowerFunc.h>
void *OPS_ElasticPowerFunc(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 5) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial ";
    opserr << "ElasticPowerFunc tag <eta> -coeff c1 c2 ... ";
    opserr << "-exp e1 e2 ... ";
    opserr << "(with at least one pair of (ci,ei) values)\n";
    return 0;
  }

  int tag[1];
  double coeffData[64];
  double expData[64];
  double eta = 0.0;
  const char *paraStr;

  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPowerFunc tag\n";
    return 0;
  }

  // check if eta is provided (odd number of inputs)
  if ((argc - 3) % 2 == 1) {
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &eta) != 0) {
      opserr << "WARNING invalid eta\n";
      opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
      return 0;
    }
    argc--;
  }

  // get coefficient values
  numData = (argc - 3) / 2;
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-coeff") == 0 || strcmp(paraStr, "-coefficient") == 0 ||
      strcmp(paraStr, "-coefficients") == 0) {
    if (OPS_GetDoubleInput(&numData, coeffData) != 0) {
      opserr << "WARNING invalid coefficients\n";
      opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -coeff but got " << paraStr << endln;
    opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
    return 0;
  }
  Vector coefficients(coeffData, numData);

  // get exponent values
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-exp") == 0 || strcmp(paraStr, "-exponent") == 0 ||
      strcmp(paraStr, "-exponents") == 0) {
    if (OPS_GetDoubleInput(&numData, expData) != 0) {
      opserr << "WARNING invalid exponents\n";
      opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -exp but got " << paraStr << endln;
    opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
    return 0;
  }
  Vector exponents(expData, numData);

  // Parsing was successful, allocate the material
  theMaterial = new ElasticPowerFunc(tag[0], coefficients, exponents, eta);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ";
    opserr << "ElasticPowerFunc\n";
    return 0;
  }

  return theMaterial;
}
