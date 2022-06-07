

#include <g3_api.h>


#include <SRC/material/uniaxial/SmoothPSConcrete.h>
void *OPS_SmoothPSConcrete(G3_Runtime *rt)
{
  int argc = OPS_GetNumRemainingInputArgs() + 2;
  if (argc < 6 || argc > 9) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: uniaxialMaterial SmoothPSConcrete tag? fc? fu? Ec? "
              "<eps0?> <epsu?> <eta?>\n";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete tag\n";
    return 0;
  }

  // double fu, Ec, fc;
  double data[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr
        << "WARNING invalid uniaxialMaterial SmoothPSConcrete double inputs\n";
    return 0;
  }

  // double eps0=0.002;
  // double epsu=0.005;
  // double eta=0.2;
  double opt[3] = {0.002, 0.005, 0.2};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 3)
    numdata = 3;
  if (OPS_GetDoubleInput(&numdata, opt) < 0) {
    opserr
        << "WARNING invalid uniaxialMaterial SmoothPSConcrete double inputs\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  return new SmoothPSConcrete(tag, data[0], data[1], data[2], opt[0], opt[1],
                              opt[2]);
}
