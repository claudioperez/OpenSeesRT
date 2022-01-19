

#include <g3_api.h>


#include <SRC/element/XMUelements/ASI3D8QuadWithSensitivity.h>
void *OPS_ASID8QuadWithSensitivity(void)
{

  static int idData[9];

  // if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 9) {
    opserr << "element ASI3D8Quad incorrect num args .. 9 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element ASI3D8Quad error reading first 9 integers\n";
    return 0;
  }

  Element *theEle = new ASI3D8QuadWithSensitivity(
      idData[0], idData[1], idData[2], idData[3], idData[4], idData[5],
      idData[6], idData[7], idData[8]);
  return theEle;
}
