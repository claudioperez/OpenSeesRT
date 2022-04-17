

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthImpact3D.h>
void *OPS_ZeroLengthImpact3D(void)
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyZeroLengthImpact3D == 0) {
    opserr << "Using ZeroLengthImpact3D element - Developed by Prof. Arash E. "
              "Zaghi & Majid Cashany, University of Connecticut (UConn) "
              "Copyright 2012 - Use at your Own Peril\n";
    numMyZeroLengthImpact3D++;
  }

  // get the id and end nodes
  int iData[4];
  double dData[7];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D tag\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D 1st node " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D 2nd node " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D direction " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D initial gap input " << eleTag
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D frictionRatio " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Ktangent " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Knormal " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[4]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Kn2 Input " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[5]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D Delta_y Input " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D cohesion " << eleTag << endln;
    return 0;
  }

  /*
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D origin X " << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[8]) != 0) {
    opserr << "WARNING ZeroLengthImpact3D origin Y " << eleTag << endln;
    return 0;
  }
  */

  // now create the ZeroLengthImpact3D and add it to the Domain

  Element *theZeroLengthImpact3D = new ZeroLengthImpact3D(
      eleTag, iData[1], iData[2], iData[3], dData[0], dData[1], dData[2],
      dData[3], dData[4], dData[5], dData[6]);

  if (theZeroLengthImpact3D == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    // delete theMaterial;
    return 0;
  }

  return theZeroLengthImpact3D;
}
