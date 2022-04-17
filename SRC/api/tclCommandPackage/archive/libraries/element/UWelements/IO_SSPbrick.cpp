

#include <g3_api.h>


#include <SRC/element/UWelements/SSPbrick.h>
OPS_Export void *OPS_SSPbrick(void)
{
  if (num_SSPbrick == 0) {
    num_SSPbrick++;
    opserr << "SSPbrick element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 10) {
    opserr
        << "Invalid #args, want: element SSPbrick eleTag? iNode? jNode? kNode? "
           "lNode? mNode? nNode? pNode? qNode? matTag? <b1? b2? b3?>\n";
    return 0;
  }

  int iData[10];
  double dData[3];
  dData[0] = 0.0;
  dData[1] = 0.0;
  dData[2] = 0.0;

  int numData = 10;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPbrick " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[9];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SSPbrick " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  if (numRemainingInputArgs == 13) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid optional data: element SSPbrick " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPbrick(iData[0], iData[1], iData[2], iData[3], iData[4],
                            iData[5], iData[6], iData[7], iData[8],
                            *theMaterial, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPbrick\n";
    return 0;
  }

  return theElement;
}
