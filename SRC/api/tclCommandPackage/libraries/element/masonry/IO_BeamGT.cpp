

#include <g3_api.h>


#include <SRC/element/masonry/BeamGT.h>
void *OPS_BeamGT()

{
  // print out a message about who wrote this element & any copyright info
  // wanted
  // if (numMyBeam == 0) {
  //	opserr << " \n";
  //    opserr << "           Beam with Flexure and Shear Hinges\n";
  //    opserr << "   Written by Gonzalo Torrisi UNCuyo Copyright 2016\n";
  //	opserr << "                 Only in plane X-Y \n";
  //    opserr << "                Use at your Own Peril\n";

  //  numMyBeam++;
  //}

  Element *theBeam = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theBeam = new BeamGT();
    return theBeam;
  }

  if (numRemainingArgs != 14) {
    opserr << "ERROR - BeamGT not enough args provided, want: element BeamGT "
              "tag? Node1? Node2?  matTag? matTag2? matTag3? E? G? A? I? Lp1? "
              "Lp2? Lr? fc?\n";
    // numMyBeam++;
  }

  // get the id and end nodes
  int iData[6];
  double dData[8];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING error reading element material 1 tag for element "
           << eleTag << endln;
    return 0;
  }
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING error reading element material 2 tag for element "
           << eleTag << endln;
    return 0;
  }
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING error reading element material 3 tag for element "
           << eleTag << endln;
    return 0;
  }
  int matID = iData[3];
  int matID2 = iData[4];
  int matID3 = iData[5];

  numData = 8;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading Elastic properties for element" << eleTag
           << endln;
    return 0;
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
  UniaxialMaterial *theMaterial2 = OPS_GetUniaxialMaterial(matID2);
  UniaxialMaterial *theMaterial3 = OPS_GetUniaxialMaterial(matID3);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element "
           << eleTag << endln;
    return 0;
  }

  if (theMaterial2 == 0) {
    opserr << "WARNING material with tag " << matID2 << "not found for element "
           << eleTag << endln;
    return 0;
  }
  if (theMaterial3 == 0) {
    opserr << "WARNING material with tag " << matID3 << "not found for element "
           << eleTag << endln;
    return 0;
  }
  // now create the truss and add it to the Domain

  theBeam = new BeamGT(eleTag, iData[1], iData[2], *theMaterial, *theMaterial2,
                       *theMaterial3, dData[0], dData[1], dData[2], dData[3],
                       dData[4], dData[5], dData[6], dData[7]);

  if (theBeam == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    delete theMaterial;
    delete theMaterial2;
    delete theMaterial3;
    return 0;
  }

  return theBeam;
}
