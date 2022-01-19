

#include <g3_api.h>


#include <SRC/material/uniaxial/PenaltyMaterial.h>
OPS_Export void *OPS_PenaltyMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING insufficient args, uniaxialMaterial Penalty $tag "
              "$otherTag $penalty"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Penalty $tag $otherTag $penalty"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial Penalty tag: "
           << iData[0] << endln;
    return 0;
  }

  double penalty = 0.0;
  numData = 1;
  if (OPS_GetDouble(&numData, &penalty) < 0) {
    opserr << "WARNING invalid input uniaxialMaterial Penalty tag: " << iData[0]
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new PenaltyMaterial(iData[0], *theOtherMaterial, penalty);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "PenaltyMaterial\n";
    return 0;
  }

  return theMaterial;
}
