

#include <g3_api.h>


#include <SRC/element/dispBeamColumn/DispBeamColumn3dWithSensitivity.h>
void *OPS_DispBeamColumn3dWithSensitivity()
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr
        << "insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag "
           "<-mass mass> <-cmass>\n";
    return 0;
  }

  // inputs:
  int iData[5];
  int numData = 5;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  // options
  double mass = 0.0;
  numData = 1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0) {
          opserr << "WARNING: invalid mass\n";
          return 0;
        }
      }
    }
  }

  // check transf
  CrdTransf *theTransf = OPS_getCrdTransf(iData[3]);
  if (theTransf == 0) {
    opserr << "coord transfomration not found\n";
    return 0;
  }

  // check beam integrataion
  BeamIntegrationRule *theRule = OPS_getBeamIntegrationRule(iData[4]);
  if (theRule == 0) {
    opserr << "beam integration not found\n";
    return 0;
  }
  BeamIntegration *bi = theRule->getBeamIntegration();
  if (bi == 0) {
    opserr << "beam integration is null\n";
    return 0;
  }

  // check sections
  const ID &secTags = theRule->getSectionTags();
  SectionForceDeformation **sections =
      new SectionForceDeformation *[secTags.Size()];
  for (int i = 0; i < secTags.Size(); i++) {
    sections[i] = OPS_getSectionForceDeformation(secTags(i));
    if (sections[i] == 0) {
      opserr << "section " << secTags(i) << "not found\n";
      delete[] sections;
      return 0;
    }
  }

  Element *theEle = new DispBeamColumn3dWithSensitivity(
      iData[0], iData[1], iData[2], secTags.Size(), sections, *bi, *theTransf,
      mass);
  delete[] sections;
  return theEle;
}
