

#include <g3_api.h>


#include <SRC/element/UWelements/Brick8FiberOverlay.h>
void *OPS_Brick8FiberOverlay(void)
{
  if (num_Brick8FiberOverlay == 0) {
    num_Brick8FiberOverlay++;
    opserr << "Brick8FiberOverlay element - Written: M.Chiaramonte, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 15) {
    opserr << "Want: Brick8FiberOverlay tag? nd1? nd2? nd3? nd4? nd5? nd6? "
              "nd7? nd8? matTag? AreaFiber? B1? B2? B3? B4?\n";
    return 0;
  }

  int iData[9];
  double dData[5];
  int matTag = 0;
  int eleTag = 0;
  int numData = 9;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Brick8FiberOverlay"
           << endln;
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Brick8FiberOverlay: invalid matTag for element: "
           << eleTag << "\n";
    return 0;
  }
  numData = 5;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Brick8FiberOverlay " << eleTag
           << endln;
    return 0;
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element "
           << eleTag << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new Brick8FiberOverlay(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], iData[5], iData[6], iData[7],
                                      iData[8], *theMaterial, dData[0],
                                      dData[1], dData[2], dData[3], dData[4]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Brick8FiberOverlay\n";
    return 0;
  }

  return theElement;
}
