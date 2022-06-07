

#include <g3_api.h>


#include <SRC/element/UWelements/SSPbrickUP.h>
OPS_Export void *OPS_SSPbrickUP(void)
{
  if (num_SSPbrickUP == 0) {
    num_SSPbrickUP++;
    opserr << "SSPbrickUP element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 17) {
    opserr << "Invalid #args, want: element SSPbrickUP eleTag? iNode? jNode? "
              "kNode? lNode? mNode? nNode? pNode? qNode? matTag? fBulk? fDen? "
              "k1? k2? k3? e? alpha? <b1? b2? b3?>\n";
    return 0;
  }

  int iData[10];
  double dData[10];
  dData[7] = 0.0;
  dData[8] = 0.0;
  dData[9] = 0.0;

  int numData = 10;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPbrickUP " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[9];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SSPbrickUP " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element SSPbrickUP " << iData[0]
           << endln;
    return 0;
  }

  if (numRemainingInputArgs == 20) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPbrickUP " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPbrickUP(
      iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6],
      iData[7], iData[8], *theMaterial, dData[0], dData[1], dData[2], dData[3],
      dData[4], dData[5], dData[6], dData[7], dData[8], dData[9]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPbrickUP\n";
    return 0;
  }

  return theElement;
}
