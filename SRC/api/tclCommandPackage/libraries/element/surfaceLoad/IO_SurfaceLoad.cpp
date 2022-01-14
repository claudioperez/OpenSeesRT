

#include <g3_api.h>


#include <SRC/element/surfaceLoad/SurfaceLoad.h>
void *OPS_SurfaceLoad(void)
{
  if (num_SurfaceLoad == 0) {
    num_SurfaceLoad++;
    opserr << "SurfaceLoad element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 6) {
    opserr << "Want: element SurfaceLoad eleTag?  iNode? jNode? kNode? lNode? "
              "pressure?\n";
    return 0;
  }

  int iData[5];
  double dData[1];

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SurfaceLoadElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SurfaceLoad " << iData[0] << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new SurfaceLoad(iData[0], iData[1], iData[2], iData[3], iData[4],
                               dData[0]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SurfaceLoadElement\n";
    return 0;
  }

  return theElement;
}
