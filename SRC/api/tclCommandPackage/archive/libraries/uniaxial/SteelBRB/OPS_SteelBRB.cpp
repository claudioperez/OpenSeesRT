

#include <g3_api.h>


#include <SRC/material/uniaxial/SteelBRB.h>
void *OPS_SteelBRB(void)
{

  int tag;
  double E, sigmaY0, sigmaY_T, alpha_T, beta_T, delta_T, sigmaY_C, alpha_C,
      beta_C, delta_C, Tol;

  Tol = 1.0e-14;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 11 && numArgs != 12) { // ---- correct!!
    opserr << "Warning Insufficient args: unixialMaterial SteelBRB tag E "
              "sigmaY0 sigmaY_T alpha_T beta_T delta_T sigmaY_C alpha_C beta_C "
              "delta_C <Tol> \n";
    return 0;
  }

  int iData[1];
  double dData[11];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial SimplifiedJ2 \n";
    return 0;
  }
  tag = iData[0];

  numData = numArgs - 1;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid double values: nDMaterial SimplifiedJ2 " << tag
           << endln;
    return 0;
  }

  E = dData[0];
  sigmaY0 = dData[1];
  sigmaY_T = dData[2];
  alpha_T = dData[3];
  beta_T = dData[4];
  delta_T = dData[5];
  sigmaY_C = dData[6];
  alpha_C = dData[7];
  beta_C = dData[8];
  delta_C = dData[9];
  if (numArgs == 12) {
    Tol = dData[10];
  }

  // Parsing was successful, allocate the material
  UniaxialMaterial *theMaterial =
      new SteelBRB(tag, E, sigmaY0, sigmaY_T, alpha_T, alpha_C, sigmaY_C,
                   beta_T, beta_C, delta_T, delta_C, Tol);
  return theMaterial;
}
