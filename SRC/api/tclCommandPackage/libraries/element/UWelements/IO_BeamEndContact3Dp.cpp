

#include <g3_api.h>


#include <SRC/element/UWelements/BeamEndContact3Dp.h>
OPS_Export void *OPS_BeamEndContact3Dp(void)
{
  if (num_BeamEndContact3Dp == 0) {
    num_BeamEndContact3Dp++;
    opserr << "BeamEndContact3Dp element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 6) {
    opserr << "Invalid #args, want: element BeamEndContact3Dp eleTag? iNode? "
              "jNode? sNode? radius? penalty? <cFlag>?\n";
    return 0;
  }

  int iData[4];
  double dData[2];
  int icSwitch = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamEndContact3Dp "
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element BeamEndContact3Dp "
           << iData[0] << endln;
    return 0;
  }

  numRemainingInputArgs -= 6;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr
          << "WARNING invalid initial contact flag: element BeamEndContact3Dp "
          << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamEndContact3Dp(iData[0], iData[1], iData[2], iData[3],
                                     dData[0], dData[1], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type "
              "BeamEndContact3DpElement\n";
    return 0;
  }

  return theElement;
}
