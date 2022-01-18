

#include <g3_api.h>


#include <SRC/material/uniaxial/ViscousMaterial.h>
void *OPS_ViscousMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3 || numArgs > 4) {
    opserr << "Invalid #args,  want: uniaxialMaterial Viscous tag? C? alpha? "
              "<minVel?> ... "
           << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  dData[2] = 1.0e-11; // setting default minVel

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Viscous" << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Viscous " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new ViscousMaterial(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Viscous\n";
    return 0;
  }

  return theMaterial;
}
