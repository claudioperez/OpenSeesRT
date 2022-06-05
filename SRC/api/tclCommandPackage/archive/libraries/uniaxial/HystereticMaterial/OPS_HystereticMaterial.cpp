

#include <g3_api.h>


#include <SRC/material/uniaxial/HystereticMaterial.h>
void *OPS_HystereticMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 18 && numArgs != 17 && numArgs != 14 && numArgs != 13) {
    opserr << "Want: uniaxialMaterial Hysteretic tag? mom1p? rot1p? mom2p? "
              "rot2p? <mom3p? rot3p?> "
           << "\nmom1n? rot1n? mom2n? rot2n? <mom3n? rot3n?> pinchX? pinchY? "
              "damfc1? damfc2? <beta?>";
    return 0;
  }

  int iData[1];
  double dData[17];
  for (int i = 0; i < 17; i++)
    dData[i] = 0.0;

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Hysteretic" << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Hysteretic " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  if (numData > 13)
    theMaterial = new HystereticMaterial(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16]);
  else
    theMaterial =
        new HystereticMaterial(iData[0], dData[0], dData[1], dData[2], dData[3],
                               dData[4], dData[5], dData[6], dData[7], dData[8],
                               dData[9], dData[10], dData[11], dData[12]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Hysteretic\n";
    return 0;
  }

  return theMaterial;
}
