

#include <g3_api.h>


#include <SRC/material/uniaxial/SteelECThermal.h>
void *OPS_SteelECThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[2];
  double dData[6];
  const char *typeChar = new char[20];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelECThermal tag?" << endln;
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() == 3 ||
      OPS_GetNumRemainingInputArgs() == 7) {
    typeChar = OPS_GetString();

    if (strcmp(typeChar, "EC3") == 0) {
      iData[1] = 3;
    } else if (strcmp(typeChar, "EC2Nh") == 0 ||
               strcmp(typeChar, "EC2NH") == 0) {
      iData[1] = 21;
    } else if (strcmp(typeChar, "EC2NC") == 0 ||
               strcmp(typeChar, "EC2Nc") == 0) {
      iData[1] = 22;
    } else if (strcmp(typeChar, "EC2X") == 0 || strcmp(typeChar, "EC2x") == 0) {
      iData[1] = 23;
    } else {
      opserr << "WARNING invalid material type for uniaxialMaterial "
                "SteelECThermal "
             << iData[0] << endln;
      return 0;
    }

  } else if (OPS_GetNumRemainingInputArgs() == 2 ||
             OPS_GetNumRemainingInputArgs() == 6) {
    iData[1] = 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 2 && numData != 6) {
    opserr << "Invalid #args, want: uniaxialMaterial SteelECThermal "
           << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial SteelECThermal "
           << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 2) {
    dData[2] = STEEL_01_DEFAULT_A1;
    dData[3] = STEEL_01_DEFAULT_A2;
    dData[4] = STEEL_01_DEFAULT_A3;
    dData[5] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SteelECThermal(iData[0], iData[1], dData[0], dData[1],
                                   dData[2], dData[3], dData[4], dData[5]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "SteelECThermal Material\n";
    return 0;
  }

  return theMaterial;
}
