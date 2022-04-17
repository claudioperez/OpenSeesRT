

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthND.h>
void *OPS_ZeroLengthND()
{
  int ndm = OPS_GetNDM();
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING too few arguments "
           << "want - element zeroLengthND eleTag? iNode? jNode? "
           << "NDTag? <1DTag?>"
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";

    return 0;
  }

  int idata[4];
  numdata = 4;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: failed to get integer data\n";
    return 0;
  }
  NDMaterial *nmat = OPS_getNDMaterial(idata[3]);
  if (nmat == 0) {
    opserr << "WARNING: NDMaterial " << idata[3] << " is not defined\n";
    return 0;
  }

  UniaxialMaterial *umat = 0;
  int uniTag;
  if (OPS_GetIntInput(&numdata, &uniTag) >= 0) {
    umat = OPS_getUniaxialMaterial(uniTag);
    if (umat == 0) {
      opserr << "WARNING: uniaxial material " << uniTag << " is not defined\n";
      return 0;
    }
  } else {
    OPS_ResetCurrentInputArg(-1);
  }

  const char *type = OPS_GetString();
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;
  if (strcmp(type, "-orient") == 0) {
    if (OPS_GetNumRemainingInputArgs() < 6) {
      opserr << "WARNING: insufficient orient values\n";
      return 0;
    }
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
      opserr << "WARNING: invalid double input\n";
      return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
      opserr << "WARNING: invalid double input\n";
      return 0;
    }
  }

  if (umat == 0) {
    return new ZeroLengthND(idata[0], ndm, idata[1], idata[2], x, y, *nmat);
  } else {
    return new ZeroLengthND(idata[0], ndm, idata[1], idata[2], x, y, *nmat,
                            *umat);
  }
}
