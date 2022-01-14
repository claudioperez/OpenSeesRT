

#include <g3_api.h>


#include <SRC/element/UWelements/SSPquad.h>
OPS_Export void *OPS_SSPquad(void)
{
  if (num_SSPquad == 0) {
    num_SSPquad++;
    opserr << "SSPquad element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args, want: element SSPquad eleTag? iNode? jNode? "
              "kNode? lNode? matTag? type? thickness? <b1? b2?>?\n";
    return 0;
  }

  int iData[6];
  const char *theType;
  double dData[3] = {1.0, 0.0, 0.0};

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPquad " << iData[0]
           << endln;
    return 0;
  }

  theType = OPS_GetString();

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid thickness data: element SSPquad " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SSPquad " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  if (numRemainingInputArgs == 10) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquad " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPquad(iData[0], iData[1], iData[2], iData[3], iData[4],
                           *theMaterial, theType, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPquad\n";
    return 0;
  }

  return theElement;
}
