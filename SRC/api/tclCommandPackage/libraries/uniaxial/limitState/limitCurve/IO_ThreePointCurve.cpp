

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/limitCurve/ThreePointCurve.h>
void *OPS_ThreePointCurve(G3_Runtime *rt)
{
  if (OPS_GetNumRemainingInputArgs() < 12) {
    opserr << "WARNING insufficient arguments\n";
    opserr
        << "Want: limitCurve ThreePoint tag? eleTag? x1? y1? x2? y2? x3? y3?";
    opserr << "Kdeg? Fres? defType? forType?" << endln;
    opserr << "<ndI? ndJ? dof? perpDirn?>" << endln;
    return 0;
  }
  int tag;
  int eleTag;
  double Kdeg;
  double Fres;
  int defType, forType;
  double x1, y1;
  double x2, y2;
  double x3, y3;
  int ndI = 0;
  int ndJ = 0;
  int dof = 0;
  int perpDirn = 0;

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid limitCurve ThreePoint tag" << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
    opserr << "WARNING invalid element tag for associated beam-column element "
              "(eleTag)\n";
    opserr << "LimitCurve ThreePoint: " << tag << endln;
    return 0;
  }

  numdata = 8;
  double data[8];
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double data\n";
    opserr << "limitCurve ThreePoint: " << tag << endln;
    return 0;
  }
  x1 = data[0];
  y1 = data[1];
  x2 = data[2];
  y2 = data[3];
  x3 = data[4];
  y3 = data[5];
  Kdeg = data[6];
  Fres = data[7];

  numdata = 1;
  if (OPS_GetIntInput(&numdata, &defType) < 0) {
    opserr << "WARNING invalid deformation type defType\n";
    opserr << "LimitCurve ThreePoint: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &forType) < 0) {
    opserr << "WARNING invalid force type forType\n";
    opserr << "LimitCurve ThreePoint: " << tag << endln;
    return 0;
  }
  if (defType == 2) {

    if (OPS_GetNumRemainingInputArgs() < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr
          << "Want: limitCurve ThreePoint tag? eleTag? x1? y1? x2? y2? x3? y3?";
      opserr << "Kdeg? Fres? defType? forType?" << endln;
      opserr << "ndI? ndJ? dof? perpDirn?" << endln;
    }
    if (OPS_GetIntInput(&numdata, &ndI) < 0) {
      opserr << "WARNING invalid node I\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &ndJ) < 0) {
      opserr << "WARNING invalid node J\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &dof) < 0) {
      opserr << "WARNING invalid degree of freedom for drift\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
      opserr << "WARNING invalid direction for column length\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return 0;
    }
  }
  // Parsing was successful, allocate the material
  // Subtract one from dof and perpDirn for C indexing
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;
  return new ThreePointCurve(tag, eleTag, theDomain, x1, y1, x2, y2, x3, y3,
                             Kdeg, Fres, defType, forType, ndI, ndJ, dof - 1,
                             perpDirn - 1);
}
