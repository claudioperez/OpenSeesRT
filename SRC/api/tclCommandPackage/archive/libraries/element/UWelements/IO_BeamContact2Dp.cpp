

#include <g3_api.h>


#include <SRC/element/UWelements/BeamContact2Dp.h>
OPS_Export void *OPS_BeamContact2Dp(void)
{
  if (num_BeamContact2Dp == 0) {
    num_BeamContact2Dp++;
    opserr << "BeamContact2Dp element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args, want: element BeamContact2Dp eleTag? iNode? "
              "jNode? secondaryNode? matTag? width? penalty? <cSwitch>?\n";
    return 0;
  }

  int iData[5];
  double dData[2];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact2Dp "
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact2Dp " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[4];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact2Dp " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 7;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact2Dp "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact2Dp(iData[0], iData[1], iData[2], iData[3],
                                  *theMaterial, dData[0], dData[1], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact2Dp\n";
    return 0;
  }

  return theElement;
}
