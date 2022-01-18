

#include <g3_api.h>


#include <SRC/element/zeroLength/CoupledZeroLength.h>
void *OPS_CoupledZeroLength()
{

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new CoupledZeroLength();
    return theEle;
  }

  if (numRemainingArgs != 6 && numRemainingArgs != 7) {
    opserr << "ERROR - CoupledZeroLength not enough args provided, want: "
              "element CoupledZeroLength tag? iNode? jNode? dirn1? dirn2? "
              "matTag? <useRayleigh?>\n";
  }

  // get the id and end nodes
  int iData[7];
  int numData;

  iData[6] = 0; // turn off rayleigh damping by default

  numData = numRemainingArgs;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];
  int matID = iData[5];
  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element "
           << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain

  theEle = new CoupledZeroLength(eleTag, iData[1], iData[2], *theMaterial,
                                 iData[3] - 1, iData[4] - 1, iData[6]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    delete theMaterial;
    return 0;
  }

  return theEle;
}
