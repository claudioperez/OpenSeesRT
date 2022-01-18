

#include <g3_api.h>


#include <SRC/element/UWelements/BeamEndContact3D.h>
OPS_Export void *OPS_BeamEndContact3D(void)
{
  if (num_BeamEndContact3D == 0) {
    num_BeamEndContact3D++;
    opserr << "BeamEndContact3D element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args, want: element BeamEndContact3D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? radius? gapTol? forceTol "
              "<cFlag>?\n";
    return 0;
  }

  int iData[5];
  double dData[3];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamEndContact3D "
           << iData[0] << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element BeamEndContact3D "
           << iData[0] << endln;
    return 0;
  }

  numRemainingInputArgs -= 8;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr
          << "WARNING invalid initial contact flag: element BeamEndContact3D "
          << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamEndContact3D(iData[0], iData[1], iData[2], iData[3], iData[4],
                           dData[0], dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type BeamEndContact3DElement\n";
    return 0;
  }

  return theElement;
}
