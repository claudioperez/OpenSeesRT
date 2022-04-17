

#include <g3_api.h>


#include <SRC/element/UWelements/BeamContact3D.h>
OPS_Export void *OPS_BeamContact3D(void)
{
  if (num_BeamContact3D == 0) {
    num_BeamContact3D++;
    // OPS_Error("BeamContact3D element - Written: K.Petek, P.Arduino,
    // P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "BeamContact3D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 10) {
    opserr << "Invalid #args,  want: element BeamContact3D eleTag?  iNode? "
              "jNode? secondaryNode? lambdaNode? radius? crdTransf? matTag? "
              "tolGap? tolF? <cSwitch>?\n";
    return 0;
  }

  int iData[7];
  double dData[3];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement"
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0]
           << endln;
    return 0;
  }

  int transfTag = iData[5];
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[6];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 10;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact3D "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact3D(iData[0], iData[1], iData[2], iData[3],
                                 iData[4], dData[0], *theTransf, *theMaterial,
                                 dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact3DElement\n";
    return 0;
  }

  return theElement;
}
