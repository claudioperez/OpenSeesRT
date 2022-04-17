

#include <g3_api.h>


#include <SRC/element/XMUelements/VS3D4QuadWithSensitivity.h>
void *OPS_VS3D4WuadWithSensitivity(void)
{

  double _rho = 1.0;
  double _R = 1.0;
  double _alphaN = 1.33;
  double _alphaT = 0.67;

  static int idData[5];
  static double dData[6];

  dData[2] = _rho;
  dData[3] = _R;
  dData[4] = _alphaN;
  dData[5] = _alphaT;

  // if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 9 || argc > 11) {
    opserr << "element Vs3D4 incorrect num args .. between 9 and 11 expected\n";
    return 0;
  }

  int numData = 5;
  if (OPS_GetIntInput(&numData, idData) != 0) {
    opserr << "element Vs3D4 error reading first 5 integers\n";
    return 0;
  }

  numData = argc - 5;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "element Vs3D4 error reading last few doubles for element"
           << idData[0] << endln;
    return 0;
  }

  return new VS3D4QuadWithSensitivity(idData[0], idData[1], idData[2],
                                      idData[3], idData[4], dData[0], dData[1],
                                      dData[2], dData[3], dData[4], dData[5]);
}
