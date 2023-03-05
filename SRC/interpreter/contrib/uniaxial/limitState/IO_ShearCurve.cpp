

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/limitCurve/ShearCurve.h>
void *OPS_ShearCurve(G3_Runtime *rt)
{
  if (OPS_GetNumRemainingInputArgs() < 12) {
    opserr << "WARNING insufficient arguments\n";
    //	    printCommand(argc,argv); // Commented out by Terje
    opserr
        << "Want: limitCurve Shear tag? eleTag? rho? fc? b? h? d? Fsw? "; // SDK
    opserr << "Kdeg? Fres? defType? forType?" << endln;
    opserr << "<ndI? ndJ? dof? perpDirn? delta?>" << endln;
    return 0;
  }
  int tag;
  int eleTag;
  double Kdeg;
  double Fres;
  int defType, forType;
  double rho;
  double fc;
  double b, h, d;
  int ndI = 0;
  int ndJ = 0;
  int dof = 0;
  int perpDirn = 0;
  double Fsw = 0.0; // SDK
  double delta = 0.0;

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid limitCurve Shear tag" << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
    opserr << "WARNING invalid element tag for associated beam-column element "
              "(eleTag)\n";
    opserr << "LimitCurve Shear: " << tag << endln;
    return 0;
  }

  numdata = 8;
  double data[8];
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    opserr << "limitCurve Shear: " << tag << endln;
    return 0;
  }
  rho = data[0];
  fc = data[1];
  b = data[2];
  h = data[3];
  d = data[4];
  Fsw = data[5];
  Kdeg = data[6];
  Fres = data[7];

  numdata = 1;
  if (OPS_GetIntInput(&numdata, &defType) < 0) {
    opserr << "WARNING invalid deformation type defType\n";
    opserr << "LimitCurve Shear: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &forType) < 0) {
    opserr << "WARNING invalid force type forType\n";
    opserr << "LimitCurve Shear: " << tag << endln;
    return 0;
  }
  if (defType == 2) {

    if (OPS_GetNumRemainingInputArgs() < 4) {
      opserr << "WARNING insufficient arguments\n";
      //	    printCommand(argc,argv); // Commented out by Terje
      opserr
          << "Want: limitCurve Shear tag? eleTag? rho? fc? b? h? d? Fsw? "; // SDK
      opserr << "Kdeg? Fres? defType? forType?" << endln;
      opserr << "ndI? ndJ? dof? perpDirn? <delta?>" << endln;
      return 0;
    }

    if (OPS_GetIntInput(&numdata, &ndI) < 0) {
      opserr << "WARNING invalid node I\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &ndJ) < 0) {
      opserr << "WARNING invalid node J\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &dof) < 0) {
      opserr << "WARNING invalid degree of freedom for drift\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
      opserr << "WARNING invalid direction for column length\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return 0;
    }
  }

  if (OPS_GetNumRemainingInputArgs() > 0) {
    if (OPS_GetDoubleInput(&numdata, &delta) < 0) {
      opserr << "WARNING invalid shift in drift surface (delta)\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return 0;
    }
  }

  // Parsing was successful, allocate the material
  // Subtract one from dof and perpDirn for C indexing
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;
  return new ShearCurve(tag, eleTag, theDomain, rho, fc, b, h, d, Fsw, Kdeg,
                        Fres, defType, forType, // SDK
                        ndI, ndJ, dof - 1, perpDirn - 1, delta);
}
