

#include <g3_api.h>


#include <SRC/material/uniaxial/Elastic2Material.h>
void *OPS_Elastic2(G3_Runtime *rt)
{
  int argc = OPS_GetNumRemainingInputArgs() + 2;

  if (argc < 4 || argc > 5) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>\n";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial Elastic tag\n";
    return 0;
  }

  // E, eta
  double data[2] = {0.0, 0.0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 2)
    numdata = 2;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  return new Elastic2Material(tag, data[0], data[1]);
}
