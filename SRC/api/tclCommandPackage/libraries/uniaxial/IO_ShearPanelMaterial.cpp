

#include <g3_api.h>


#include <SRC/material/uniaxial/ShearPanelMaterial.h>
void *OPS_ShearPanelMaterial(G3_Runtime *rt)
{
  int argc = OPS_GetNumRemainingInputArgs() + 2;
  if (argc != 42 && argc != 31) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial ShearPanel tag? stress1p? strain1p? "
              "stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
           << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? strain3n? "
              "stress4n? strain4n?> rDispP? rForceP? uForceP? "
           << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? "
              "gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
           << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? gammaFLimit? "
              "gammaE? YieldStress? ";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid uniaxialMaterial ShearPanel tag\n";
    return 0;
  }

  // double stress1p, stress2p, stress3p, stress4p;
  // double strain1p, strain2p, strain3p, strain4p;
  double datap[8];
  numdata = 8;
  if (OPS_GetDoubleInput(&numdata, datap) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  // double stress1n, stress2n, stress3n, stress4n;
  // double strain1n, strain2n, strain3n, strain4n;
  double datan[8];
  numdata = 8;
  if (argc == 42) {
    if (OPS_GetDoubleInput(&numdata, datan) < 0) {
      opserr << "WARNING invalid double inputs\n";
      return 0;
    }
  }

  // double rDispP, rForceP, uForceP
  double dataP[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, dataP) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  // rDispN, rForceN, uForceN;
  double dataN[3];
  numdata = 3;
  if (argc == 42) {
    if (OPS_GetDoubleInput(&numdata, dataN) < 0) {
      opserr << "WARNING invalid double inputs\n";
      return 0;
    }
  }

  // double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
  // double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
  // double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
  // double gammaE, yStr;
  double data[17];
  numdata = 17;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  // allocate the pinching material
  if (argc == 42) {
    return new ShearPanelMaterial(
        tag, datap[0], datap[1], datap[2], datap[3], datap[4], datap[5],
        datap[6], datap[7], datan[0], datan[1], datan[2], datan[3], datan[4],
        datan[5], datan[6], datan[7], dataP[0], dataP[1], dataP[2], dataN[0],
        dataN[1], dataN[2], data[0], data[0], data[0], data[0], data[0],
        data[0], data[0], data[0], data[0], data[0], data[0], data[0], data[0],
        data[0], data[0], data[15], data[16]);
  }
  if (argc == 31) {
    return new ShearPanelMaterial(
        tag, datap[0], datap[1], datap[2], datap[3], datap[4], datap[5],
        datap[6], datap[7], dataP[0], dataP[1], dataP[2], data[0], data[0],
        data[0], data[0], data[0], data[0], data[0], data[0], data[0], data[0],
        data[0], data[0], data[0], data[0], data[0], data[15], data[16]);
  }
}
