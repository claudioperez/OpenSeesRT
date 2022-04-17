

#include <g3_api.h>


#include <SRC/element/UWelements/BeamContact2D.h>
OPS_Export void *OPS_BeamContact2D(void)
{
  if (num_BeamContact2D == 0) {
    num_BeamContact2D++;
    opserr << "BeamContact2D element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 9) {
    opserr << "Invalid #args, want: element BeamContact2D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? matTag? width? gapTol? "
              "forceTol? <cSwitch>?\n";
    return 0;
  }

  int iData[6];
  double dData[3];
  int icSwitch = 0;

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact2D " << iData[0]
           << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact2D " << dData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact2D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 9;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact2D "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamContact2D(iData[0], iData[1], iData[2], iData[3], iData[4],
                        *theMaterial, dData[0], dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact2DElement\n";
    return 0;
  }

  return theElement;
}
