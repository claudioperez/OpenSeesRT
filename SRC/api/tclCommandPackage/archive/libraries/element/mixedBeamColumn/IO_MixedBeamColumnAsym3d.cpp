

#include <g3_api.h>


#include <SRC/element/mixedBeamColumn/MixedBeamColumnAsym3d.h>
void *OPS_MixedBeamColumnAsym3d()
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
  int cmass = 0;
  double dData[2]; // input of ys and zs
  dData[0] = 0.0;
  dData[1] = 0.0;
  int doRayleigh = 1;
  bool geomLinear = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-cMass") == 0) {
      opserr << "WARNING: consistent mass not implemented\n";
    } else if (strcmp(type, "-mass") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0) {
          opserr << "WARNING: invalid mass\n";
          return 0;
        }
      }
    } else if (strcmp(type, "-shearCenter") == 0) {
      // Get the coordinates of shear center w.r.t centroid
      numData = 2;
      if (OPS_GetDoubleInput(&numData, dData) < 0) {
        opserr << "WARNING: invalid ys and zs\n";
        return 0;
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
        opserr
            << "WARNING: Invalid doRayleigh in element MixedBeamColumnAsym3d "
            << iData[0];
        return 0;
      }
    } else if (strcmp(type, "-geomLinear") == 0) {
      opserr
          << "WARNING: geometric linear in the basic system not implemented\n";
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

  Element *theEle = new MixedBeamColumnAsym3d(
      iData[0], iData[1], iData[2], secTags.Size(), sections, *bi, *theTransf,
      dData[0], dData[1], mass, doRayleigh, geomLinear);
  delete[] sections;
  return theEle;
}
void *OPS_MixedBeamColumnAsym3dTcl()
{
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  double dData2[2]; // input of ys and zs
  dData2[0] = 0.0;
  dData2[1] = 0.0;
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != 3) {
    opserr << "ERROR: MixedBeamColumnAsym3d: invalid number of dimensions\n";
    return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != 6) {
    opserr << "ERROR: MixedBeamColumnAsym3d: invalid number of degrees of "
              "freedom\n";
    return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "ERROR: MixedBeamColumnAsym3d: too few arguments\n";
    return 0;
  }

  // Get required input data
  numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - MixedBeamColumnAsym3d\n";
    return 0;
  }
  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int numIntgrPts = iData[3];
  int secTag = iData[4];
  int transfTag = iData[5];

  // Get the section
  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(secTag);
  if (theSection == 0) {
    opserr << "WARNING section with tag " << secTag << "not found for element "
           << eleTag << endln;
    return 0;
  }

  SectionForceDeformation **sections =
      new SectionForceDeformation *[numIntgrPts];
  for (int i = 0; i < numIntgrPts; i++) {
    sections[i] = theSection;
  }

  // Get the coordinate transformation
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag
           << "not found for element " << eleTag << endln;
    return 0;
  }

  // Set Default Values for Optional Input
  int doRayleigh = 1;
  double massDens = 0.0;
  bool geomLinear = false;
  BeamIntegration *beamIntegr = 0;

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

    } else if (strcmp(sData, "-integration") == 0) {
      const char *sData2 = OPS_GetString();

      if (strcmp(sData2, "Lobatto") == 0) {
        beamIntegr = new LobattoBeamIntegration();
      } else if (strcmp(sData2, "Legendre") == 0) {
        beamIntegr = new LegendreBeamIntegration();
      } else if (strcmp(sData2, "Radau") == 0) {
        beamIntegr = new RadauBeamIntegration();
      } else if (strcmp(sData2, "NewtonCotes") == 0) {
        beamIntegr = new NewtonCotesBeamIntegration();
      } else if (strcmp(sData2, "Trapezoidal") == 0) {
        beamIntegr = new TrapezoidalBeamIntegration();
      } else if (strcmp(sData2, "RegularizedLobatto") == 0 ||
                 strcmp(sData2, "RegLobatto") == 0) {
        numData = 4;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING invalid input, want: -integration "
                    "RegularizedLobatto $lpI $lpJ $zetaI $zetaJ \n";
          return 0;
        }
        BeamIntegration *otherBeamInt = 0;
        otherBeamInt = new LobattoBeamIntegration();
        beamIntegr = new RegularizedHingeIntegration(
            *otherBeamInt, dData[0], dData[1], dData[2], dData[3]);
        if (otherBeamInt != 0) {
          delete otherBeamInt;
        }
      } else {
        opserr << "WARNING invalid integration type, element: " << eleTag;
        return 0;
      }
    } else if (strcmp(sData, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
        opserr
            << "WARNING: Invalid doRayleigh in element MixedBeamColumnAsym3d "
            << eleTag;
        return 0;
      }

    } else if (strcmp(sData, "-geomLinear") == 0) {
      geomLinear = true;

    } else if (strcmp(sData, "-shearCenter") == 0) {
      // Get the coordinates of shear center w.r.t centroid
      numData = 2;
      if (OPS_GetDoubleInput(&numData, &dData2[0]) < 0) {
        opserr << "WARNING: invalid ys and zs\n";
        return 0;
      }
    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  // Set the beam integration object if not in options
  if (beamIntegr == 0) {
    beamIntegr = new LobattoBeamIntegration();
  }

  // now create the element and add it to the Domain
  Element *theElement = new MixedBeamColumnAsym3d(
      eleTag, nodeI, nodeJ, numIntgrPts, sections, *beamIntegr, *theTransf,
      dData2[0], dData2[1], massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  delete[] sections;
  if (beamIntegr != 0)
    delete beamIntegr;

  return theElement;
}
