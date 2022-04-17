

#include <g3_api.h>


#include <SRC/element/UWelements/PileToe3D.h>
OPS_Export void *OPS_PileToe3D(void)
{
  if (num_PileToe3D == 0) {
    num_PileToe3D++;
    // OPS_Error("PileToe3D element - Written: P.Arduino, P.Mackenzie-Helnwein,
    // U.Washington\n", 1);
    opserr << "PileToe3D element - Written: P.Arduino, P.Mackenzie-Helnwein, "
              "U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args,  want: element PileToe3D eleTag?  iNode? BiNode? "
              "BjNode? radius? k? crdTransf?\n";
    return 0;
  }

  int iData[5];
  double dData[2];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element PileToe3D" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid radius data: element PileToe3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid  k data: element PileToe3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer crdTransf data: element PileToe3D"
           << iData[0] << endln;
    return 0;
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element PileToe3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the element
  theElement = new PileToe3D(iData[0], iData[1], iData[2], iData[3], dData[0],
                             dData[1], *theTransf);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type PileToe3D\n";
    return 0;
  }

  return theElement;
}
