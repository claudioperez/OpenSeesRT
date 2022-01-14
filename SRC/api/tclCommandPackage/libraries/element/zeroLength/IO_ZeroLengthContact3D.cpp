

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthContact3D.h>
void *OPS_ZeroLengthContact3D()
{
  // a quick check on number of args
  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "ZeroLengthContact3D::WARNING too few arguments "
           << "want - element ZeroLengthContact3D eleTag? iNode? jNode? Kn? "
              "Kt? fs? c? dir?";
    return 0;
  }

  // get the ele tag
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "ZeroLengthContact3D::WARNING invalied int inputs\n";
    return 0;
  }
  int eleTag = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];

  // read the material properties
  double ddata[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "ZeroLengthContact3D::WARNING invalied double inputs\n";
    return 0;
  }
  double Kn = ddata[0];
  double Kt = ddata[1];
  double fs = ddata[2];
  double c = ddata[3];

  int dir;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &dir) < 0) {
    opserr << "ZeroLengthContact3D::WARNING invalied direction\n";
    return 0;
  }

  //
  // now we create the element and add it to the domain
  //

  double originX, originY;
  originX = 0;
  originY = 0;

  if (dir == 0) {
    if (OPS_GetNumRemainingInputArgs() > 1) {
      if (OPS_GetDoubleInput(&numdata, &originX) < 0) {
        opserr << "ZeroLengthContact3D::WARNING invalied originX\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&numdata, &originY) < 0) {
        opserr << "ZeroLengthContact3D::WARNING invalied originY\n";
        return 0;
      }
    }
  }

  return new ZeroLengthContact3D(eleTag, iNode, jNode, dir, Kn, Kt, fs, c,
                                 originX, originY);
}
