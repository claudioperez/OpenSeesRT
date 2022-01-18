

#include <g3_api.h>


#include <SRC/material/uniaxial/UniaxialJ2Plasticity.h>
void *OPS_UniaxialJ2Plasticity(G3_Runtime *rt)
{
  int argc = OPS_GetNumRemainingInputArgs() + 2;
  if (argc < 7) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: uniaxialMaterial UniaxialJ2Plasticity tag? E? sigmaY? "
              "Hkin? <Hiso?>\n";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial UniaxialJ2Plasticity tag\n";
    return 0;
  }

  // double E, sigmaY, Hkin, Hiso;
  // Hiso =0.0;
  double data[4] = {0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 4)
    numdata = 4;

  // Parsing was successful, allocate the material
  return new UniaxialJ2Plasticity(tag, data[0], data[1], data[2], data[3]);
}
