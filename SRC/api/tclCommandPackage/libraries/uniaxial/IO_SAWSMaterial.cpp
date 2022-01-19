

#include <g3_api.h>


#include <SRC/material/uniaxial/SAWSMaterial.h>
OPS_Export void *OPS_SAWSMaterial(void)
{
  if (numSAWSMaterials == 0) {
    numSAWSMaterials++;
    // OPS_Error("SAWSMaterial unaxial material - Written by Paxti Uriz,
    // Exponent 2009\n", 1);
    opserr << "SAWSMaterial unaxial material - Written by Paxti Uriz, Exponent "
              "2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[10];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SAWSMaterial tag" << endln;
    return 0;
  }

  numData = 10;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial SAWS tag? F0? FI? dU? S0?"
           << endln;
    opserr << "    R1? R2? R3? R4? alpha? beta?" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SAWSMaterial(iData[0], dData[0], dData[1], dData[2],
                                 dData[3], dData[4], dData[5], dData[6],
                                 dData[7], dData[8], dData[9]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type SAWSMaterial\n";
    return 0;
  }

  return theMaterial;
}
