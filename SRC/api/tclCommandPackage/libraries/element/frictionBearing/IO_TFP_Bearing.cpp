

#include <g3_api.h>


#include <SRC/element/frictionBearing/TFP_Bearing.h>
void *OPS_TFP_Bearing()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyBearing == 0) {
    opserr << "TFP_Bearing element - Written by Tracy Becker, UC Berkeley "
              "Copyright 2011\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new TFP_Bearing();
    return theEle;
  }

  if (numRemainingArgs != 25 && numRemainingArgs != 24 &&
      numRemainingArgs != 26 != numRemainingArgs != 27) {
    opserr << "ERROR - TFP_Bearing incorrect # args provided, want: element "
              "TFP_Bearing tag? iNode? jNode? ";
    opserr << "$R1 $R2 $R3 $R4 $do1 $do2 $do3 $do4 $din1 $din2 $din3 $din4 "
              "$mu1 $mu2 $mu3 $mu4";
    opserr << " $h1 $h2 $h3 $h4 $H0 <$a> <$K>\n";
    return theEle;
  }

  // get the id and end nodes
  int iData[3];
  double dData[24];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  if (numRemainingArgs == 24) {
    numData = 21;
    dData[21] = 10.0; // initial Axial Load = 0.0
    dData[22] = 1.0e12;
    dData[23] = 0.01;
  } else if (numRemainingArgs == 25) {
    numData = 22;
    dData[22] = 1.0e12;
    dData[23] = 0.01;
  } else if (numRemainingArgs == 26) {
    numData = 23;
    dData[22] = 1.0e12;
  } else {
    numData = 24;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element area for element" << eleTag
           << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  int ndm = OPS_GetNDM();
  if (ndm == 3) {
    theEle = new TFP_Bearing(eleTag, iData[1], iData[2], &dData[0], &dData[4],
                             &dData[8], &dData[12], &dData[16], dData[20],
                             dData[21], dData[23], dData[22]);
  } else {
    theEle = new TFP_Bearing2d(eleTag, iData[1], iData[2], &dData[0], &dData[4],
                               &dData[8], &dData[12], &dData[16], dData[20],
                               dData[21], dData[23], dData[22]);
  }

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  return theEle;
}
