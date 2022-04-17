

#include <g3_api.h>


#include <SRC/element/frictionBearing/TPB1D.h>
void *OPS_TPB1D()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyTPB1D == 0) {
    opserr << "TPB1D2D element - Written by Troy/Fenz UC Berkeley Copyright "
              "2011 - Use at your Own Peril\n";
    numMyTPB1D++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new TPB1D();
    return theEle;
  }

  if (numRemainingArgs != 20) {
    opserr << "ERROR - TPB1D2D not enough args provided, want: element TPB1D2D "
              "tag? iNode? jNode? direction? mu1? mu2? mu3? R1? R2? R3? h1? "
              "h2? h3? D1? D2? D3? d1? d2? d3? W?\n";
    numMyTPB1D++;
  }

  // get the id and end nodes
  int iData[4];
  double dData[16];
  int numData;

  numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 16;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element area for element" << eleTag
           << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  theEle = new TPB1D(eleTag, iData[1], iData[2], iData[3] - 1, &dData[0],
                     &dData[3], &dData[6], &dData[9], &dData[12], dData[15]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element type TPB1D with tag "
           << eleTag << endln;
    return 0;
  }

  return theEle;
}
