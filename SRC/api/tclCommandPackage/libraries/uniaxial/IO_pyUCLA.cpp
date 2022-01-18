

#include <g3_api.h>


#include <SRC/material/uniaxial/pyUCLA.h>
OPS_Export void *OPS_pyUCLA(void)
{
  if (num_pyUCLA == 0) {
    num_pyUCLA++;
    // OPS_Error("pyUCLAMaterial unaxial material - Written by H.Shin,
    // P.Arduino, U.Washington\n based on model of E.Taciroglu, C.Rha,
    // J.Wallace, UCLA", 1);
    opserr << "pyUCLAMaterial unaxial material - Written by H.Shin, P.Arduino, "
              "U.Washington\n based on model of E.Taciroglu, C.Rha, J.Wallace, "
              "UCLA\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() != 5) {
    opserr << "Invalid #args,  want: uniaxialMaterial pyUCLA tag? soilType? "
              "pult? y50? Cd? "
           << endln;
    return 0;
  }

  int iData[2];
  double dData[3];

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag or soilType uniaxialMaterial pyUCLAMaterial"
           << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid pyData data for material uniaxial pyUCLA " << iData[0]
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new pyUCLA(iData[0], iData[1], dData[0], dData[1], dData[2]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type pyUCLAMaterial\n";
    return 0;
  }

  return theMaterial;
}
