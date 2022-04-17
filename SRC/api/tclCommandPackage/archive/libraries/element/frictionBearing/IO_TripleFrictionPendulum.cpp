

#include <g3_api.h>


#include <SRC/element/frictionBearing/TripleFrictionPendulum.h>
void *OPS_TripleFrictionPendulum()
{
  if (numTripleFrictionPendulum == 0) {
    numTripleFrictionPendulum++;
    opserr << "TripleFrictionPendulum element v2.0.0 - Written by Nhan@unr\n";
  }

  // get the id and end nodes
  int iData[10];
  double dData[11];
  int numData, eleTag;

  numData = 10;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data";
    return 0;
  }
  eleTag = iData[0];

  // get the friction models
  FrictionModel *theFrnMdls[3];
  for (int i = 0; i < 3; i++) {
    theFrnMdls[i] = OPS_getFrictionModel(iData[3 + i]);
    if (theFrnMdls[i] == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << iData[3 + i] << endln;
      opserr << "TripleFrictionPendulum element: " << eleTag << endln;
      return 0;
    }
  }

  // get the uniaxial materials
  UniaxialMaterial *theMaterials[4];
  for (int i = 0; i < 4; i++) {
    theMaterials[i] = OPS_getUniaxialMaterial(iData[6 + i]);
    if (theMaterials[i] == 0) {
      opserr << "WARNING uniaxial material not found\n";
      opserr << "uniaxialMaterial: " << iData[6 + i] << endln;
      opserr << "TripleFrictionPendulum element: " << eleTag << endln;
      return 0;
    }
  }

  numData = 11;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element" << eleTag << endln;
    return 0;
  }

  // create the element and add it to the Domain
  Element *theTripleFrictionPendulum = new TripleFrictionPendulum(
      eleTag, iData[1], iData[2], theFrnMdls, theMaterials, dData[0], dData[1],
      dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8],
      dData[9], dData[10]);

  if (theTripleFrictionPendulum == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  return theTripleFrictionPendulum;
}
