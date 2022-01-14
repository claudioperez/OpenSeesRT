

#include <g3_api.h>


#include <SRC/element/UWelements/BeamContact3Dp.h>
OPS_Export void *OPS_BeamContact3Dp(void)
{
  if (num_BeamContact3Dp == 0) {
    num_BeamContact3Dp++;
    // OPS_Error("BeamContact3Dp element - Written: K.Petek, C.McGann,
    // P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "BeamContact3Dp element - Written: K.Petek, C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args,  want: element BeamContact3Dp eleTag?  iNode? "
              "jNode? secondaryNode? radius? crdTransf? matTag? penalty? "
              "<cSwitch>?\n";
    return 0;
  }

  int iData[6];
  double dData[2];
  int icSwitch = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0]
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement"
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0]
           << endln;
    return 0;
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 8;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact3Dp "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamContact3Dp(iData[0], iData[1], iData[2], iData[3], dData[0],
                         *theTransf, *theMaterial, dData[1], icSwitch);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type BeamContact3DpElement\n";
    return 0;
  }

  return theElement;
}
