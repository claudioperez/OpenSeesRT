

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticMultiLinear.h>
void *OPS_ElasticMultiLinear(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 7) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial ";
    opserr << "ElasticMultiLinear tag <eta> -strain strainPoints ";
    opserr << "-stress stressPoints  ";
    opserr << "(with at least two stress-strain points)\n";
    return 0;
  }

  int tag[1];
  double strainData[64];
  double stressData[64];
  double eta = 0.0;
  const char *paraStr;

  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticMultiLinear tag\n";
    return 0;
  }

  // check if eta is provided (odd number of inputs)
  if ((argc - 3) % 2 == 1) {
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &eta) != 0) {
      opserr << "WARNING invalid eta\n";
      opserr << "uniaxialMaterial ElasticMultiLinear: " << tag[0] << endln;
      return 0;
    }
    argc--;
  }

  // get strain data points
  numData = (argc - 3) / 2;
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-strain") == 0) {
    if (OPS_GetDoubleInput(&numData, strainData) != 0) {
      opserr << "WARNING invalid strainPoints\n";
      opserr << "uniaxialMaterial ElasticMultiLinear: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -strain but got " << paraStr << endln;
    opserr << "uniaxialMaterial ElasticMultiLinear: " << tag[0] << endln;
    return 0;
  }
  Vector strainPts(strainData, numData);

  // get stress data points
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-stress") == 0) {
    if (OPS_GetDoubleInput(&numData, stressData) != 0) {
      opserr << "WARNING invalid stressPoints\n";
      opserr << "uniaxialMaterial ElasticMultiLinear: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -stress but got " << paraStr << endln;
    opserr << "uniaxialMaterial ElasticMultiLinear: " << tag[0] << endln;
    return 0;
  }
  Vector stressPts(stressData, numData);

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMultiLinear(tag[0], strainPts, stressPts, eta);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ";
    opserr << "ElasticMultiLinear\n";
    return 0;
  }

  return theMaterial;
}
