

#include <g3_api.h>


#include <SRC/element/UWelements/SSPquadUP.h>
OPS_Export void *OPS_SSPquadUP(void)
{
  if (num_SSPquadUP == 0) {
    num_SSPquadUP++;
    opserr << "SSPquadUP element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();
  // LM change
  if (numRemainingInputArgs < 13) {
    opserr << "Invalid #args, want: element SSPquadUP eleTag? iNode? jNode? "
              "kNode? lNode? matTag? t? fBulk? fDen? k1? k2? e? alpha? <b1? "
              "b2?> <Pup? Plow? Pleft? Pright?>?\n";
    return 0;
  }

  int iData[6];
  double dData[13];
  dData[7] = 0.0;
  dData[8] = 0.0;
  dData[9] = 0.0;
  dData[10] = 0.0;
  dData[11] = 0.0;
  dData[12] = 0.0;
  // LM change

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPquadUP " << iData[0]
           << endln;
    return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element SSPquadUP " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SSPquadUP " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // LM change
  if (numRemainingInputArgs == 15) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0]
             << endln;
      return 0;
    }
  } else if (numRemainingInputArgs == 19) {
    numData = 6;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPquadUP(
      iData[0], iData[1], iData[2], iData[3], iData[4], *theMaterial, dData[0],
      dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
      dData[8], dData[9], dData[10], dData[11], dData[12]);
  // LM change

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPquadUP\n";
    return 0;
  }

  return theElement;
}
