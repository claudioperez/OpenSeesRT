#include <g3_api.h>
#include <Vector.h>
#include <SRC/element/zeroLength/ZeroLengthRocking.h>

void *OPS_ZeroLengthRocking()
{

  int ndm = OPS_GetNDM(); // the spatial dimension of the problem

  //
  // first scan the command line to obtain eleID, iNode, jNode, and the
  // orientation of ele xPrime and yPrime not along the global x and y axis
  //

  // a quick check on number of args
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING too few arguments "
           << "want - element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";

    return 0;
  }

  // eleTag, iNode, jNode
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalied int inputs "
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return 0;
  }
  int eleTag = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];

  // look for rocking required inputs
  double ddata[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalied double inputs "
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return 0;
  }
  double kr = ddata[0];
  double R = ddata[1];
  double theta = ddata[2];
  double kap = ddata[3];

  // create the vectors for the element orientation
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;
  double xi = 1.0e-8;
  double dTol = 1.0e-7;
  double vTol = 1.0e-7;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *flag = OPS_GetString();
    numdata = 1;
    if (strcmp(flag, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr
            << "WARNING not enough parameters after -orient flag for ele "
            << eleTag << "- element ZeroLengthRocking eleTag? iNode? jNode? "
            << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
        return 0;

      } else {
        double value;
        // read the x values
        for (int i = 0; i < 3; i++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag
                   << "- element ZeroLength eleTag? iNode? jNode? "
                   << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? "
                      "y3?>\n";
            return 0;
          } else {
            x(i) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag
                   << "- element ZeroLength eleTag? iNode? jNode? "
                   << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? "
                      "y3?>\n";
            return 0;
          } else {
            y(j) = value;
          }
        }
      }

    } else if (strcmp(flag, "-xi") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING not enough parameters after -xi flag for ele "
               << eleTag << endln;
        return 0;
      } else {
        if (OPS_GetDoubleInput(&numdata, &xi) < 0) {
          opserr << "WARNING invalid -xi value for ele  " << eleTag << endln;
          return 0;
        }
      }

    } else if (strcmp(flag, "-dTol") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING not enough parameters after -dTol flag for ele "
               << eleTag << endln;
        return 0;
      } else {
        if (OPS_GetDoubleInput(&numdata, &dTol) < 0) {
          opserr << "WARNING invalid -dTol value for ele  " << eleTag << endln;
          return 0;
        }
      }

    } else if (strcmp(flag, "-vTol") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING not enough parameters after -vTol flag for ele "
               << eleTag << endln;
        return 0;
      } else {
        if (OPS_GetDoubleInput(&numdata, &vTol) < 0) {
          opserr << "WARNING invalid -vTol value for ele  " << eleTag << endln;
          return 0;
        }
      }
    }
  }

  //
  // now we create the element and add it to the domain
  //
  return new ZeroLengthRocking(eleTag, ndm, iNode, jNode, x, y, kr, R, theta,
                               kap, xi, dTol, vTol);
}
