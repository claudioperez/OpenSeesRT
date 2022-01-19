

#include <g3_api.h>


#include <SRC/element/UWelements/SimpleContact2D.h>
OPS_Export void *OPS_SimpleContact2D(void)
{
  if (num_SimpleContact2D == 0) {
    num_SimpleContact2D++;
    // OPS_Error("SimpleContact2D element - Written: K.Petek, P.Arduino,
    // P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "SimpleContact2D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 8) {
    opserr << "Invalid #args,  want: element SimpleContact2D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? matTag? tolGap? tolForce?\n";
    return 0;
  }

  int iData[6];
  double dData[2];

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact2DElement"
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact2D " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact2D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new SimpleContact2D(iData[0], iData[1], iData[2], iData[3],
                                   iData[4], *theMaterial, dData[0], dData[1]);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type SimpleContact2DElement\n";
    return 0;
  }

  return theElement;
}
