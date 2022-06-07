

#include <g3_api.h>


#include <SRC/element/masonry/MasonPan12.h>
void *OPS_MasonPan12()

//#define OPS_Export

// OPS_Export void *
// OPS_MasonPan12element()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  // if (numMyPanel == 0) {
  // opserr << " \n";
  //  opserr << "                 REFINED MASONRY PANEL\n";
  //  opserr << "   Written by Gonzalo Torrisi UNCuyo Copyright 2015\n";
  //	opserr << "         Model with 6 compression struts-36 dof\n";
  //	opserr << "                 Only in plane X-Y \n";
  // opserr << "                Use at your Own Peril\n";

  //    numMyPanel++;
  //  }

  Element *thePanel = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    thePanel = new MasonPan12();
    return thePanel;
  }

  if (numRemainingArgs != 18) {
    opserr << "ERROR - Masonry Panel not enough args provided, want: element "
              "MasonryPanel tag? Node1? Node2? Node3? Node4?  Node5?  Node6?  "
              "Node7?  Node8?  Node9?   Node10?   Node11?   Node12?   matTag? "
              "matTag2? thick? wfactor? w1?\n";
    //    numMyPanel++;
  }

  // get the id and end nodes
  int iData[15];
  double dData[3];
  int numData;

  numData = 13;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[13]) != 0) {
    opserr << "WARNING error reading element material 1 tag for element "
           << eleTag << endln;
    return 0;
  }
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[14]) != 0) {
    opserr << "WARNING error reading element material 2 tag for element "
           << eleTag << endln;
    return 0;
  }

  int matID = iData[13];
  int matID2 = iData[14];

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element areas, thickness and properties "
              "for element"
           << eleTag << endln;
    return 0;
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
  UniaxialMaterial *theMaterial2 = OPS_GetUniaxialMaterial(matID2);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element "
           << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain

  thePanel = new MasonPan12(eleTag, iData[1], iData[2], iData[3], iData[4],
                            iData[5], iData[6], iData[7], iData[8], iData[9],
                            iData[10], iData[11], iData[12], *theMaterial,
                            *theMaterial2, dData[0], dData[1], dData[2]);

  if (thePanel == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    delete theMaterial;
    delete theMaterial2;
    return 0;
  }

  return thePanel;
}
