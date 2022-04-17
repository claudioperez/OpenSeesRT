

#include <g3_api.h>


#include <SRC/element/mixedBeamColumn/MixedBeamColumn3d.h>
void *OPS_MixedBeamColumn3d()
{
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  int sDataLength = 40;
  // char *sData  = new char[sDataLength];
  // char *sData2 = new char[sDataLength];
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != NDM) {
    opserr << "ERROR: MixedBeamColumn3d: invalid number of dimensions\n";
    return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
    opserr
        << "ERROR: MixedBeamColumn3d: invalid number of degrees of freedom\n";
    return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "ERROR: MixedBeamColumn3d, too few arguments: "
              "eleTag,ndI,ndJ,transfTag,integrationTag\n";
    return 0;
  }

  // Get required input data
  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - MixedBeamColumn3d\n";
    return 0;
  }
  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int transfTag = iData[3];
  int beamIntTag = iData[4];

  // Get the coordinate transformation
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag
           << "not found for element " << eleTag << endln;
    return 0;
  }

  // Get beam integrataion
  BeamIntegrationRule *theRule = OPS_getBeamIntegrationRule(beamIntTag);
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

  // Set Default Values for Optional Input
  int doRayleigh = 1;
  double massDens = 0.0;
  bool geomLinear = true;

  // Loop through remaining arguments to get optional input
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *sData = OPS_GetString();
    if (strcmp(sData, "-mass") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid input, want: -mass $massDens \n";
        return 0;
      }
      massDens = dData[0];

    } else if (strcmp(sData, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
        opserr << "WARNING: Invalid doRayleigh in element MixedBeamColumn3d "
               << eleTag;
        return 0;
      }

    } else if (strcmp(sData, "-geomNonlinear") == 0) {
      geomLinear = false;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  // now create the element and add it to the Domain
  Element *theElement =
      new MixedBeamColumn3d(eleTag, nodeI, nodeJ, secTags.Size(), sections, *bi,
                            *theTransf, massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  delete[] sections;
  return theElement;
}
