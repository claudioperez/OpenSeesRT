

#include <g3_api.h>


#include <SRC/material/uniaxial/Pinching4Material.h>
void *OPS_Pinching4Material(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata != 40 && numdata != 29) {
    opserr << "WARNING: Insufficient arguments\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    return 0;
  }

  UniaxialMaterial *mat = 0;
  int tDmg = -1;
  if (OPS_GetNumRemainingInputArgs() == 39) {
    double data[38];
    numdata = 38;
    if (OPS_GetDoubleInput(&numdata, data)) {
      return 0;
    }
    const char *type = OPS_GetString();
    if (strcmp(type, "cycle") == 0 || strcmp(type, "Cycle") == 0 ||
        strcmp(type, "DamageCycle") == 0 || strcmp(type, "damageCycle") == 0) {
      tDmg = 1;
    } else if (strcmp(type, "energy") == 0 || strcmp(type, "Energy") == 0 ||
               strcmp(type, "DamageEnergy") == 0 ||
               strcmp(type, "damageEnergy") == 0) {
      tDmg = 0;
    } else {
      opserr << "WARNING invalid type of damage calculation specified\n";
      opserr << "Pinching4 material: " << tag << endln;
      return 0;
    }

    mat = new Pinching4Material(
        tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6],
        data[7], data[8], data[9], data[10], data[11], data[12], data[13],
        data[14], data[15], data[16], data[17], data[18], data[19], data[20],
        data[21], data[22], data[23], data[24], data[25], data[26], data[27],
        data[28], data[29], data[30], data[31], data[32], data[33], data[34],
        data[35], data[36], data[37], tDmg);
  } else if (OPS_GetNumRemainingInputArgs() == 28) {
    double data[27];
    numdata = 27;
    if (OPS_GetDoubleInput(&numdata, data)) {
      return 0;
    }
    const char *type = OPS_GetString();
    if (strcmp(type, "cycle") == 0 || strcmp(type, "Cycle") == 0 ||
        strcmp(type, "DamageCycle") == 0 || strcmp(type, "damageCycle") == 0) {
      tDmg = 1;
    } else if (strcmp(type, "energy") == 0 || strcmp(type, "Energy") == 0 ||
               strcmp(type, "DamageEnergy") == 0 ||
               strcmp(type, "damageEnergy") == 0) {
      tDmg = 0;
    } else {
      opserr << "WARNING invalid type of damage calculation specified\n";
      opserr << "Pinching4 material: " << tag << endln;
      return 0;
    }

    mat = new Pinching4Material(
        tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6],
        data[7], data[8], data[9], data[10], data[11], data[12], data[13],
        data[14], data[15], data[16], data[17], data[18], data[19], data[20],
        data[21], data[22], data[23], data[24], data[25], data[26], tDmg);
  }

  if (mat == 0) {
    opserr << "WARNING: failed to create Pinching4material material\n";
    return 0;
  }

  return mat;
}
