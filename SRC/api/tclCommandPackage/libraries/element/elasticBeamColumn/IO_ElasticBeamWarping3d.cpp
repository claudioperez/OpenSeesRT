

#include <g3_api.h>


#include <SRC/element/elasticBeamColumn/ElasticBeamWarping3d.h>
void *OPS_ElasticBeamWarping3d(void)
{
  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 11 && numArgs != 6) {
    opserr << "insufficient "
              "arguments:eleTag,iNode,jNode,<A,E,G,J,Iy,Iz>or<sectionTag>,"
              "transfTag,Cw\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 7) {
    opserr << "ndm must be 3 and ndf must be 7\n";
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

  if (numArgs == 6) {
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

  // Read Cw
  numData = 1;
  double Cw;
  if (OPS_GetDoubleInput(&numData, &Cw) < 0)
    return 0;

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
    }
  }

  if (theSection != 0) {
    return new ElasticBeamWarping3d(iData[0], iData[1], iData[2], theSection,
                                    *theTrans, Cw, mass);
  } else {
    return new ElasticBeamWarping3d(iData[0], data[0], data[1], data[2],
                                    data[3], data[4], data[5], iData[1],
                                    iData[2], *theTrans, Cw, mass);
  }
}
