

#include <g3_api.h>


#include <SRC/element/surfaceLoad/TriSurfaceLoad.h>
void *OPS_TriSurfaceLoad(void)
{
  if (num_TriSurfaceLoad == 0) {
    num_TriSurfaceLoad++;
    opserr << "TriSurfaceLoad element - Written: J. A. Abell (UANDES). "
              "Inspired by the makers of SurfaceLoad\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Want: element TriSurfaceLoad eleTag?  iNode? jNode? kNode? "
              "pressure? <rhoH?>\n";
    return 0;
  }

  int iData[4];
  double dData[1];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element TriSurfaceLoadElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element TriSurfaceLoad " << iData[0]
           << endln;
    return 0;
  }

  double rhoH = 0;
  int num_args_remaining = OPS_GetNumRemainingInputArgs();

  if (num_args_remaining > 0) {
    numData = 1;
    OPS_GetDoubleInput(&numData, &rhoH);
  }

  // Parsing was successful, allocate the material
  theElement = new TriSurfaceLoad(iData[0], iData[1], iData[2], iData[3],
                                  dData[0], rhoH);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type TriSurfaceLoadElement\n";
    return 0;
  }

  return theElement;
}
