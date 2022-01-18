

#include <g3_api.h>


#include <SRC/material/uniaxial/SteelMP.h>
void *OPS_SteelMP(G3_Runtime *rt)
{
  int argc = OPS_GetNumRemainingInputArgs() + 2;

  // Check that there is the minimum number of arguments
  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial SteelMP tag? fy? E0? b? ";
    opserr << " <coeffR1?  coeffR2? a1? a2?>\n";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelMP tag\n";
    return 0;
  }

  // Read required Steel01 material parameters
  // double fy, E, b;
  double data[3];
  numdata = 3;
  if (argc < 6) {
    opserr << "WARNING insufficient number of hardening parameters\n";
    opserr << "uniaxialMaterial Steel03: " << tag << "\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelMP double inputs\n";
    return 0;
  }

  // Read optional Steel01 material parameters
  // double r, coeffR1, coeffR2, a1, a2;
  double opt[5] = {20.0, 18.5, 0.15, 0, 0};
  numdata = 5;
  if (argc > 6) {
    if (OPS_GetDoubleInput(&numdata, opt) < 0) {
      opserr << "WARNING invalid uniaxialMaterial SteelMP double inputs\n";
      return 0;
    }
  } // if

  return new SteelMP(tag, data[0], data[1], data[2], opt[0], opt[1], opt[2],
                     opt[3], opt[4]);
}
