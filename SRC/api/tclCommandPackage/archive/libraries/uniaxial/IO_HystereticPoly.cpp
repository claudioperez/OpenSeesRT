

#include <g3_api.h>


#include <SRC/material/uniaxial/HystereticPoly.h>
void *OPS_HystereticPoly(G3_Runtime *rt)
{
  UniaxialMaterial *theMaterial = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial HystereticPoly tag? ka? kb? a? b1? b2? "
              "<tol?>"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[6];
  dData[5] = 1.0e-20; // default tolerance

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial HystereticPoly"
           << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial HystereticPoly " << iData[0] << endln;
    return 0;
  }
  if (dData[0] <= 0.0) {
    opserr << "uniaxialMaterial HystereticPoly ka must be positive" << endln;
    return 0;
  }
  if (dData[1] >= dData[0]) {
    opserr << "uniaxialMaterial HystereticPoly kb must be < ka" << endln;
    return 0;
  }
  if (dData[2] <= 0.0) {
    opserr << "uniaxialMaterial HystereticPoly a must be positive and <> 1"
           << endln;
    return 0;
  }
  if (dData[2] == 1.0) {
    opserr << "uniaxialMaterial HystereticPoly a must be positive and <> 1"
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new HystereticPoly(iData[0], dData[0], dData[1], dData[2],
                                   dData[3], dData[4], dData[5]);
  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type HystereticPoly\n";
    return 0;
  }

  return theMaterial;
}
