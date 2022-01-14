

#include <g3_api.h>


#include <SRC/material/uniaxial/BarSlipMaterial.h>
void *OPS_BarSlipMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata != 15 && numdata != 13) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial BarSlip tag? ";
    opserr << "fc? fy? Es? fu? Eh? db? ld? nb? width? ";
    opserr << "depth? bsflag? type? <damage? unit?>\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    return 0;
  }

  double data[10];
  numdata = 10;
  if (OPS_GetDoubleInput(&numdata, data)) {
    return 0;
  }

  const char *bsFlag = OPS_GetString();
  int bsf = 0;
  if (strcmp(bsFlag, "strong") == 0 || strcmp(bsFlag, "Strong") == 0) {
    bsf = 0;
  } else if (strcmp(bsFlag, "weak") == 0 || strcmp(bsFlag, "Weak") == 0) {
    bsf = 1;
  } else {
    opserr << "WARNING invalid bond strength specified\n";
    opserr << "BarSlip: " << tag << "\n";
    return 0;
  }

  const char *type = OPS_GetString();
  int typ = 0;
  if (strcmp(type, "beamtop") == 0 || strcmp(type, "beamTop") == 0 ||
      strcmp(type, "beam") == 0 || strcmp(type, "Beam") == 0) {
    typ = 0;
  } else if (strcmp(type, "beambot") == 0 || strcmp(type, "beamBot") == 0 ||
             strcmp(type, "beambottom") == 0 ||
             strcmp(type, "beamBottom") == 0) {
    typ = 1;
  } else if (strcmp(type, "column") == 0 || strcmp(type, "Column") == 0) {
    typ = 2;
  } else {
    opserr << "WARNING invalid location of bar specified\n";
    opserr << "BarSlip: " << tag << "\n";
    return 0;
  }

  numdata = OPS_GetNumRemainingInputArgs();
  int dmg = -1, unt = -1;
  if (numdata > 1) {
    const char *damage = OPS_GetString();
    if (strcmp(damage, "damage1") == 0 || strcmp(damage, "Damage1") == 0) {
      dmg = 1;
    } else if (strcmp(damage, "damage2") == 0 ||
               strcmp(damage, "Damage2") == 0) {
      dmg = 2;
    } else if (strcmp(damage, "nodamage") == 0 ||
               strcmp(damage, "Nodamage") == 0 ||
               strcmp(damage, "NoDamage") == 0 ||
               strcmp(damage, "noDamage") == 0) {
      dmg = 0;
    } else {
      opserr << "WARNING invalid damage specified\n";
      opserr << "BarSlip: " << tag << "\n";
      return 0;
    }

    const char *unit = OPS_GetString();
    if (strcmp(unit, "mpa") == 0 || strcmp(unit, "MPa") == 0 ||
        strcmp(unit, "mPa") == 0 || strcmp(unit, "Mpa") == 0) {
      unt = 1;
    } else if (strcmp(unit, "psi") == 0 || strcmp(unit, "Psi") == 0 ||
               strcmp(unit, "PSI") == 0) {
      unt = 2;
    } else if (strcmp(unit, "Pa") == 0 || strcmp(unit, "pa") == 0) {
      unt = 3;
    } else if (strcmp(unit, "psf") == 0 || strcmp(unit, "Psf") == 0 ||
               strcmp(unit, "PSF") == 0) {
      unt = 4;
    } else if (strcmp(unit, "ksi") == 0 || strcmp(unit, "Ksi") == 0 ||
               strcmp(unit, "KSI") == 0) {
      unt = 5;
    } else if (strcmp(unit, "ksf") == 0 || strcmp(unit, "Ksf") == 0 ||
               strcmp(unit, "KSF") == 0) {
      unt = 6;
    } else {
      opserr << "WARNING invalid unit specified\n";
      opserr << "BarSlip: " << tag << "\n";
      return 0;
    }
  }

  UniaxialMaterial *mat = 0;
  if (dmg == -1) {
    mat = new BarSlipMaterial(tag, data[0], data[1], data[2], data[3], data[4],
                              data[5], data[6], data[7], data[8], data[9], bsf,
                              typ);
  } else {
    mat = new BarSlipMaterial(tag, data[0], data[1], data[2], data[3], data[4],
                              data[5], data[6], data[7], data[8], data[9], bsf,
                              typ, dmg, unt);
  }

  if (mat == 0) {
    opserr << "WARNING: failed to create BarSlipMaterial material\n";
    return 0;
  }

  return mat;
}
