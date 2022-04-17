

#include <g3_api.h>


#include <SRC/material/uniaxial/EPPGapMaterial.h>
void *OPS_EPPGapMaterial(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 4) {
    opserr << "Invalid #args,  want: uniaxialMaterial ElasticPPGap tag E Fy "
              "gap <eta damage>\n";
    return 0;
  }

  int tag, damage = 0;
  double dData[4];
  dData[3] = 0.0; // setting default eta to 0.

  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial EPPGap" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData > 4)
    numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial EPPGap \n";
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData > 0) {
    numData = 1;
    const char *damagestr = OPS_GetString();
    if (strcmp(damagestr, "damage") == 0 || strcmp(damagestr, "Damage") == 0) {
      damage = 1;
    }
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new EPPGapMaterial(tag, dData[0], dData[1], dData[2], dData[3], damage);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type EPPGap\n";
    return 0;
  }

  return theMaterial;
}
