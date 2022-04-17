

#include <g3_api.h>


#include <SRC/element/UWelements/SimpleContact3D.h>
OPS_Export void *OPS_SimpleContact3D(void)
{
  if (num_SimpleContact3D == 0) {
    num_SimpleContact3D++;
    // OPS_Error("SimpleContact3D element - Written: K.Petek, P.Arduino,
    // P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "SimpleContact3D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 10) {
    opserr << "Invalid #args,  want: element SimpleContact3D eleTag? iNode? "
              "jNode? kNode? lNode? secondaryNode? lambdaNode? matTag? tolGap? "
              "tolForce?\n";
    return 0;
  }

  int iData[8];
  double dData[2];

  int numData = 8;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact3DElement"
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact3D " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[7];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement =
      new SimpleContact3D(iData[0], iData[1], iData[2], iData[3], iData[4],
                          iData[5], iData[6], *theMaterial, dData[0], dData[1]);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type SimpleContact3DElement\n";
    return 0;
  }

  return theElement;
}
