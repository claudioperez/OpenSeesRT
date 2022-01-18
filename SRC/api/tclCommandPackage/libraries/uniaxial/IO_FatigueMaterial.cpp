

#include <g3_api.h>


#include <SRC/material/uniaxial/FatigueMaterial.h>
void *OPS_FatigueMaterial(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
    opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << endln;
    opserr << " <-min min?> <-max max?>" << endln;
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invlid int inputs\n";
    return 0;
  }

  double Dmax = 1.0;
  double E0 = 0.191;
  double m = -0.458;
  double epsmin = NEG_INF_STRAIN;
  double epsmax = POS_INF_STRAIN;
  numdata = 1;

  while (OPS_GetNumRemainingInputArgs() > 1) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-Dmax") == 0) {
      if (OPS_GetDouble(&numdata, &Dmax) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
      }
    } else if (strcmp(type, "-E0") == 0) {
      if (OPS_GetDouble(&numdata, &E0) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
      }
    } else if (strcmp(type, "-m") == 0) {
      if (OPS_GetDouble(&numdata, &m) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
      }
    } else if (strcmp(type, "-min") == 0) {
      if (OPS_GetDouble(&numdata, &epsmin) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
      }
    } else if (strcmp(type, "-max") == 0) {
      if (OPS_GetDouble(&numdata, &epsmax) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
      }
    }
  }

  UniaxialMaterial *mat = OPS_getUniaxialMaterial(idata[1]);
  if (mat == 0) {
    opserr << "WARNING component material does not exist\n";
    opserr << "Component material: " << idata[1];
    opserr << "\nuniaxialMaterial Fatigue: " << idata[0] << endln;
    return 0;
  }

  UniaxialMaterial *theMat =
      new FatigueMaterial(idata[0], *mat, Dmax, E0, m, epsmin, epsmax);
  if (theMat == 0) {
    opserr << "WARNING: failed to create FatigueMaterial material\n";
    return 0;
  }

  return theMat;
}
