

#include <g3_api.h>


#include <SRC/element/forceBeamColumn/ForceBeamColumn3d.h>
void *OPS_ForceBeamColumn3d()
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "insufficient "
              "arguments:eleTag,iNode,jNode,transfTag,integrationTag\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "ndm must be 3 and ndf must be 6\n";
    return 0;
  }

  // inputs:
  int iData[5];
  int numData = 5;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  // options
  double mass = 0.0, tol = 1e-12;
  int maxIter = 10;
  numData = 1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 1) {
        if (OPS_GetIntInput(&numData, &maxIter) < 0) {
          opserr << "WARNING invalid maxIter\n";
          return 0;
        }
        if (OPS_GetDoubleInput(&numData, &tol) < 0) {
          opserr << "WARNING invalid tol\n";
          return 0;
        }
      }
    } else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0) {
          opserr << "WARNING invalid mass\n";
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
      return 0;
    }
  }

  Element *theEle =
      new ForceBeamColumn3d(iData[0], iData[1], iData[2], secTags.Size(),
                            sections, *bi, *theTransf, mass, maxIter, tol);
  delete[] sections;
  return theEle;
}
