

#include <g3_api.h>


#include <SRC/material/uniaxial/StainlessECThermal.h>
void *OPS_StainlessECThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int iData[2];
  double dData[4];

  int numData = 1;
  // string gradeInput;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial StainlessECThermal tag?"
           << endln;
    return 0;
  }

  //
  // numData = OPS_GetNumRemainingInputArgs();
  // if (numData == 4 || numData == 8)
  //{
  //	//if (OPS_GetString(gradeChar, 20) == 0)
  //	//{
  //	//	opserr << "WARNING invalid gradeTag for uniaxialMaterial
  //StainlessECThermal " << iData[1] << endln;
  //	//	//return 0;
  //	// }
  //
  //	//gradeInput = OPS_GetString(gradeChar, 20);
  const char *gradeChar = OPS_GetString();
  if (strcmp(gradeChar, "Grade14301") == 0) {
    iData[1] = 1;
  } else if ((strcmp(gradeChar, "Grade14401") == 0 ||
              strcmp(gradeChar, "Grade14404") == 0)) {
    iData[1] = 2;
  } else if (strcmp(gradeChar, "Grade14571") == 0) {
    iData[1] = 3;
  } else if (strcmp(gradeChar, "Grade14003") == 0) {
    iData[1] = 4;
  } else if (strcmp(gradeChar, "Grade14462") == 0) {
    iData[1] = 5;
  } else {
    opserr << "WARNING invalid material grade for uniaxialMaterial "
              "StainlessECThermal "
           << iData[0] << endln;
    return 0;
  }
  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 3 && numData != 4) {
    opserr << "Invalid #args, want: uniaxialMaterial StainlessECThermal "
           << iData[0] << " fy? E? fu?" << endln;
    return 0;
  }
  //

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial StainlessECThermal "
           << iData[0] << " fy? E? fu?" << endln;
    return 0;
  }
  //
  if (numData == 3) {
    dData[3] = 0.0; // set Initial stress=0.0
  }

  // Parsing was successful, allocate the material
  theMaterial = new StainlessECThermal(iData[0], iData[1], dData[0], dData[1],
                                       dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "StainlessECThermal Material\n";
    return 0;
  }

  return theMaterial;
}
