

#include <g3_api.h>


#include <SRC/element/absorbentBoundaries/LysmerTriangle.h>
void *OPS_LysmerTriangle(void)
{
  if (num_LysmerTriangle == 0) {
    num_LysmerTriangle++;
    opserr << "LysmerTriangle element - Written: J. A. Abell (UANDES). "
              "www.joseabell.com\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Want: element LysmerTriangle eleTag?  iNode? jNode? kNode? rho "
              "Vp Vs? <length> <stage> \n";
    return 0;
  }

  int iData[4];
  double dData[3];
  double eleLength = 0;
  int stage = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element LysmerTriangleElement"
           << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element LysmerTriangle " << iData[0]
           << endln;
    return 0;
  }

  int num_args_remaining = OPS_GetNumRemainingInputArgs();

  // Its optional (but desireable) to input the element-length....
  if (num_args_remaining > 0) {
    numData = 1;
    OPS_GetDoubleInput(&numData, &eleLength);
  }

  // Its optional to set the element stage... will be set to 0 (damping mode) by
  // default
  if (num_args_remaining > 0) {
    numData = 1;
    OPS_GetIntInput(&numData, &stage);
  }

  // Parsing was successful, allocate the material
  theElement =
      new LysmerTriangle(iData[0], iData[1], iData[2], iData[3], dData[0],
                         dData[1], dData[2], eleLength, stage);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type LysmerTriangleElement\n";
    return 0;
  }

  return theElement;
}
