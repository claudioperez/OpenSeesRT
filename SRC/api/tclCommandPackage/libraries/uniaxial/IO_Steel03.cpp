

#include <g3_api.h>


#include <SRC/material/uniaxial/Steel03.h>
void *OPS_Steel03(G3_Runtime *rt)
{

  int argc = OPS_GetNumRemainingInputArgs() + 2;

  // Check that there is the minimum number of arguments
  if (argc < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Steel03 tag? fy? E0? b? r? cR1 cR2?";
    opserr << " <a1? a2? a3? a4?>\n";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel03 tag\n";
    return 0;
  }

  // Read required Steel01 material parameters
  // fy, E, b, r, cR1, cR2;
  double data[6];
  numdata = 6;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  // Read optional Steel01 material parameters
  // a1, a2, a3, a4
  if (argc > 9) {
    double opt[4];
    numdata = 4;

    if (argc < 13) {
      opserr << "WARNING insufficient number of hardening parameters\n";
      opserr << "uniaxialMaterial Steel03: " << tag << "\n";
      return 0;
    }
    if (OPS_GetDoubleInput(&numdata, opt) < 0) {
      opserr << "WARNING invalid double inputs\n";
      return 0;
    }

    // Parsing was successful, allocate the material
    return new Steel03(tag, data[0], data[1], data[2], data[3], data[4],
                       data[5], opt[0], opt[1], opt[2], opt[3]);
  } else
    // Parsing was successful, allocate the material
    return new Steel03(tag, data[0], data[1], data[2], data[3], data[4],
                       data[5]);
}
