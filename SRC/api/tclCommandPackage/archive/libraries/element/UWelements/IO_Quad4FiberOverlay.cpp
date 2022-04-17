

#include <g3_api.h>


#include <SRC/element/UWelements/Quad4FiberOverlay.h>
void *OPS_Quad4FiberOverlay(void)
{
  if (num_Quad4FiberOverlay == 0) {
    num_Quad4FiberOverlay++;
    opserr << "Quad4FiberOverlay element - Written: M.Chiaramonte, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 9) {
    opserr << "Want: Quad4FiberOverlay tag? nd1? nd2? nd3? nd4? matTag? "
              "CrossSectionArea? B1?  B2? \n";
    return 0;
  }

  int iData[6];
  double dData[3];
  int matTag = 0;
  int eleTag = 0;
  int numData = 0;
  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Quad4FiberOverlay"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Quad4FiberOverlay: invalid matTag for element: "
           << eleTag << "\n";
    return 0;
  }
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Quad4FiberOverlay " << eleTag
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
  theElement =
      new Quad4FiberOverlay(iData[0], iData[1], iData[2], iData[3], iData[4],
                            *theMaterial, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Quad4FiberOverlay\n";
    return 0;
  }

  return theElement;
}
