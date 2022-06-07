

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/LimitStateMaterial.h>
void *OPS_LimiStateMaterial(G3_Runtime *rt)
{
  UniaxialMaterial *mat = 0;

  int argc = OPS_GetNumRemainingInputArgs() + 2;
  if (argc != 20 && argc != 19 && argc != 16 && argc != 15 && argc != 22 &&
      argc != 23) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial LimitState tag? mom1p? rot1p? mom2p? "
              "rot2p? mom3p? rot3p? "
           << "\nmom1n? rot1n? mom2n? rot2n? mom3n? rot3n? pinchX? pinchY? "
              "damfc1? damfc2? beta? "
           << "\n<curveTag? curveType?>";
    return 0;
  }

  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double sp12[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, sp12) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double sp3[2];
  if (argc > 16) {
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, sp3) < 0) {
      opserr << "WARNING invalid double inputs\n";
      return 0;
    }
  }

  double sn12[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, sn12) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double sn3[2];
  if (argc > 16) {
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, sn3) < 0) {
      opserr << "WARNING invalid double inputs\n";
      return 0;
    }
  }

  double data[4];
  numdata = 4;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  double beta = 0.0;
  numdata = 1;
  if (argc == 20 || argc == 16 || argc >= 22) {
    if (OPS_GetDoubleInput(&numdata, &beta) < 0) {
      opserr << "WARNING invalid beta\n";
      return 0;
    }
  }

  int degrade = 0;
  if (argc == 22 || argc == 23) {

    double curveData[2];
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, curveData) < 0) {
      opserr << "WARNING invalid int inputs\n";
      return 0;
    }

    // LimitCurve *theCurve = GetLimitCurve(curveTag); //MRL Commented
    LimitCurve *theCurve = 0;                   // MRL Added
    theCurve = OPS_getLimitCurve(curveData[0]); // MRL Added

    if (theCurve == 0) {
      opserr << "WARNING limit curve does not exist\n";
      opserr << "limit curve: " << curveData[0];
      opserr << "\nLimitStateMaterial: " << tag << "\n";
      return 0;
    }

    if (argc == 23) {
      numdata = 1;
      if (OPS_GetIntInput(&numdata, &degrade) < 0) {
        opserr << "WARNING invalid degrade\n";
        return 0;
      }
    }
    mat = new LimitStateMaterial(
        tag, sp12[0], sp12[1], sp12[2], sp12[3], sp3[0], sp3[1], sn12[0],
        sn12[1], sn12[2], sn12[3], sn3[0], sn3[1], data[0], data[1], data[2],
        data[3], beta, *theCurve, curveData[1], degrade);
  }

  // Parsing was successful, allocate the material
  if (argc == 20 || argc == 19) {
    mat = new LimitStateMaterial(tag, sp12[0], sp12[1], sp12[2], sp12[3],
                                 sp3[0], sp3[1], sn12[0], sn12[1], sn12[2],
                                 sn12[3], sn3[0], sn3[1], data[0], data[1],
                                 data[2], data[3], beta);

  } else if (argc == 16 || argc == 15) {
    mat = new LimitStateMaterial(tag, sp12[0], sp12[1], sp12[2], sp12[3],
                                 sn12[0], sn12[1], sn12[2], sn12[3], data[0],
                                 data[1], data[2], data[3], beta);
  }
  return mat;
}
