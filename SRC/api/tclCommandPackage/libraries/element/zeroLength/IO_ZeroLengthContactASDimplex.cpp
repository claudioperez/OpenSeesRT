

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthContactASDimplex.h>
void *OPS_ZeroLengthContactASDimplex(void)
{

  double SmallNumber = 1.0e-6;
  Element *theElement = nullptr;

  // some kudos
  static int counter = 0;
  if (++counter == 1)
    opserr << "ZeroLengthContactASDimplex element - Implemented: Akan, OD., "
              "Petracca, M., Camata, G., Spacone, E. & Lai, CG. (2020)\n";

  // model dimension
  int ndm = OPS_GetNDM();
  if (ndm < 2 || ndm > 3) {
    opserr << "ZeroLengthContactASDimplex: Unsupported NDM (" << ndm
           << "). It should be 2 or 3\n";
    return theElement;
  }

  // a quick check on number of args
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "ZeroLengthContactASDimplex: WARNING: too few arguments \n"
           << "want - element zeroLengthContactASDimplex eleTag? iNode? jNode? "
              "Kn? Kt? mu? <-orient $x1 $x2 $x3> <-intType type?>\n";
    return theElement;
  }

  // start with mandatory inputs
  // read eleTag, iNode, jNode
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "ZeroLengthContactASDimplex: WARNING: invalid int inputs\n";
    return theElement;
  }

  // read Kn, Kt, mu
  double ddata[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "ZeroLengthContactASDimplex: WARNING: invalid double inputs\n";
    return theElement;
  }

  // continue with optional inputs
  Vector x_e(3);
  x_e(0) = 1.0;
  x_e(1) = 0.0;
  x_e(2) = 0.0;            // initialize orientation vector
  int integrationType = 0; // implicit by default
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *inputstring = OPS_GetString();
    if (strcmp(inputstring, "-orient") ==
        0) { // #1 read global element orientation
      if (ndm == 2) {
        if (OPS_GetNumRemainingInputArgs() < 2) {
          opserr << "ZeroLengthContactASDimplex: WARNING: insufficient orient "
                    "values in 2D\n";
          return theElement;
        }
        numdata = 3;
        if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
          opserr << "ZeroLengthContactASDimplex: WARNING: invalid double input "
                    "after -orient\n";
          return theElement;
        }
      } else if (ndm == 3) {
        if (OPS_GetNumRemainingInputArgs() < 3) {
          opserr << "ZeroLengthContactASDimplex: WARNING: insufficient orient "
                    "values in 3D\n";
          return theElement;
        }
        numdata = 3;
        if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
          opserr << "ZeroLengthContactASDimplex: WARNING: invalid double input "
                    "after -orient\n";
          return theElement;
        }
      } else {
        opserr << "ZeroLengthContactASDimplex: WARNING: -orient: model "
                  "dimension is invalid! \n";
        return theElement;
      }
    } else if (strcmp(inputstring, "-intType") ==
               0) { // #2 read type of integration
      numdata = 1;
      if (OPS_GetIntInput(&numdata, &integrationType) < 0) {
        opserr << "ZeroLengthContactASDimplex: WARNING: invalid integer after "
                  "-intType\n";
        return theElement;
      }
    }
  }
  // input reading stage is complete

  // check integration type and pick implicit if neither 1 or 0
  if (integrationType != 1 && integrationType != 0) {
    opserr << "ZeroLengthContactASDimplex: WARNING: type of integration is set "
              "to IMPLICIT due to invalid flag\n";
    integrationType = false;
  }
  // check the normal vector and normalize
  if (x_e.Norm() < SmallNumber) {
    opserr << "ZeroLengthContactASDimplex: WARNING: normal vector is NOT "
              "valid!: -orient $x1 $x2 $x3 cannot be < 0, 0, 0 >\n";
    return theElement;
  }
  x_e.Normalize(); // normalized it on input!

  // finally, create the element
  theElement = new ZeroLengthContactASDimplex(
      idata[0], idata[1], idata[2], ddata[0], ddata[1], ddata[2], ndm,
      integrationType, x_e[0], x_e[1], x_e[2]);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element zeroLengthContactASDimplex "
           << idata[0]
           << " iNode? jNode? Kn? Kt? mu? <-orient $x1 $x2 $x3> <-intType "
              "type?>\n";
  }

  return theElement;
}
