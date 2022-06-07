

#include <g3_api.h>


#include <SRC/element/elasticBeamColumn/ElasticBeam3d.h>
void *OPS_ElasticBeam3d(void)
{
  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 10 && numArgs != 5) {
    opserr << "insufficient "
              "arguments:eleTag,iNode,jNode,<A,E,G,J,Iy,Iz>or<sectionTag>,"
              "transfTag\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "ndm must be 3 and ndf must be 6\n";
    return 0;
  }

  // inputs:
  int iData[3];
  int numData = 3;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0)
    return 0;

  SectionForceDeformation *theSection = 0;
  CrdTransf *theTrans = 0;
  double data[6];
  int transfTag, secTag;
  int releasez = 0;
  int releasey = 0;

  if (numArgs == 5) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &secTag) < 0)
      return 0;
    if (OPS_GetIntInput(&numData, &transfTag) < 0)
      return 0;

    theSection = OPS_getSectionForceDeformation(secTag);
    if (theSection == 0) {
      opserr << "no section is found\n";
      return 0;
    }
    theTrans = OPS_getCrdTransf(transfTag);
    if (theTrans == 0) {
      opserr << "no CrdTransf is found\n";
      return 0;
    }
  } else {
    numData = 6;
    if (OPS_GetDoubleInput(&numData, &data[0]) < 0)
      return 0;
    numData = 1;
    if (OPS_GetIntInput(&numData, &transfTag) < 0)
      return 0;
    theTrans = OPS_getCrdTransf(transfTag);
    if (theTrans == 0) {
      opserr << "no CrdTransf is found\n";
      return 0;
    }
  }

  // options
  double mass = 0.0;
  int cMass = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string theType = OPS_GetString();
    if (theType == "-mass") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0)
          return 0;
      }
    } else if (theType == "-cMass") {
      cMass = 1;
    } else if (theType == "-releasez") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &releasez) < 0) {
          opserr << "WARNING: failed to get releasez";
          return 0;
        }
      }
    } else if (theType == "-releasey") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &releasey) < 0) {
          opserr << "WARNING: failed to get releasey";
          return 0;
        }
      }
    }
  }

  if (theSection != 0) {
    return new ElasticBeam3d(iData[0], iData[1], iData[2], theSection,
                             *theTrans, mass, cMass, releasez, releasey);
  } else {
    return new ElasticBeam3d(iData[0], data[0], data[1], data[2], data[3],
                             data[4], data[5], iData[1], iData[2], *theTrans,
                             mass, cMass, releasez, releasey);
  }
}
