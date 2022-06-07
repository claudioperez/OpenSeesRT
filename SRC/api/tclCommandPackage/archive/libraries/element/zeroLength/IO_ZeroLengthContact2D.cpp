

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthContact2D.h>
void *OPS_ZeroLengthContact2D()
{
  if (OPS_GetNumRemainingInputArgs() < 9) {
    opserr << "ZeroLengthContact2D::WARNING too few arguments "
           << "want - element ZeroLengthContact2D eleTag? iNode? jNode? Kn? "
              "Kt? fs? -normal Nx? Ny?";
    return 0;
  }

  // eleTag, iNode, jNode;
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  // Kn, Kt, fs
  double data[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING: invalid double inputs\n";
    return 0;
  }

  const char *type = OPS_GetString();
  if (strcmp(type, "-normal") != 0) {
    opserr << "ZeroLengthContact2D:: expecting "
           << "- element ZeroLengthContact2D eleTag? iNode? jNode? Kn? Kt? fs? "
              "-normal Nx? Ny? \n";
    return 0;
  }

  Vector normaldir(2);
  numdata = 2;
  if (OPS_GetDoubleInput(&numdata, &normaldir(0)) < 0) {
    opserr << "WARNING: invalid double inputs\n";
    return 0;
  }

  return new ZeroLengthContact2D(idata[0], idata[1], idata[2], data[0], data[1],
                                 data[2], normaldir);
}
