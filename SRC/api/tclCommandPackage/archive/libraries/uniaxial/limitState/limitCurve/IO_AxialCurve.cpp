

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/limitCurve/AxialCurve.h>
void *OPS_AxialCurve(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 7) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: limitCurve Axial tag? eleTag? Fsw? Kdeg? Fres? defType? "
              "forType?"
           << endln; // SDK
    opserr << "<ndI? ndJ? dof? perpDirn? delta? eleRemove?>" << endln;
    return 0;
  }
  int tag;
  int eleTag;
  double Fsw; // SDK
  double Kdeg;
  double Fres;
  int defType, forType;
  int ndI = 0;
  int ndJ = 0;
  int dof = 0;
  int perpDirn = 0;
  int eleRemove = 0;
  double delta = 0.0;

  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid Axial LimitCurve tag" << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
    opserr << "WARNING invalid element tag for associated beam-column element "
              "(eleTag)\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (OPS_GetDoubleInput(&numdata, &Fsw) < 0) {
    opserr << "WARNING invalid Fsw\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (OPS_GetDoubleInput(&numdata, &Kdeg) < 0) {
    opserr << "WARNING invalid degrading slope Kdeg\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (OPS_GetDoubleInput(&numdata, &Fres) < 0) {
    opserr << "WARNING invalid residual capacity Fres\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &defType) < 0) {
    opserr << "WARNING invalid deformation type defType\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &forType) < 0) {
    opserr << "WARNING invalid force type forType\n";
    opserr << "LimitCurve Axial: " << tag << endln;
    return 0;
  }
  if (defType == 2) {
    if (OPS_GetNumRemainingInputArgs() < 4) {
      opserr << "WARNING invalid number of arguments\n";
      opserr << "Want: limitCurve Axial tag? eleTag? Fsw? Kdeg? Fres? defType? "
                "forType?"
             << endln; // SDK
      opserr << "ndI? ndJ? dof? perpDirn? <delta? eleRemove?>" << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &ndI) < 0) {
      opserr << "WARNING invalid node I\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &ndJ) < 0) {
      opserr << "WARNING invalid node J\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &dof) < 0) {
      opserr << "WARNING invalid degree of freedom for drift\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
      opserr << "WARNING invalid direction for column length\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
  }
  if (OPS_GetNumRemainingInputArgs() > 0) {
    if (OPS_GetDoubleInput(&numdata, &delta) < 0) {
      opserr << "WARNING invalid shift in drift surface (delta)\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
  }
  if (OPS_GetNumRemainingInputArgs() > 0) {
    if (OPS_GetIntInput(&numdata, &eleRemove) < 0) {
      opserr << "WARNING invalid element removal option\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return 0;
    }
  }
  // Parsing was successful, allocate the limit curve
  // Subtract one from dof and perpDirn for C indexing
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;
  // return new AxialCurve(interp, tag, eleTag, theDomain, Fsw, //SDK
  //  			  Kdeg, Fres, defType, forType, ndI, ndJ, dof-1,
  //  perpDirn-1, 			  delta, eleRemove);
}
