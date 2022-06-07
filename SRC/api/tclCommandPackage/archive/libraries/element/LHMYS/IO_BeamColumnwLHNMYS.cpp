

#include <g3_api.h>


#include <SRC/element/LHMYS/BeamColumnwLHNMYS.h>
void *OPS_BeamColumnwLHNMYS(void)
{
  if (OPS_GetNumRemainingInputArgs() < 11) {
    opserr << "insufficient arguments:eleTag,iNode,jNode,A,E,Iz,NpI NpJ,MpI "
              "MpJ,transfTag\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 2 || ndf != 3) {
    opserr << "ndm must be 2 and ndf must be 3\n";
    return 0;
  }

  // inputs:
  int iData[3];
  int numData = 3;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING failed to read integers\n";
    return 0;
  }

  double data[7];
  // Read A, E, I, NpI, NpJ, MpI, MpJ
  numData = 7;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING failed to read doubles\n";
    return 0;
  }

  numData = 1;
  int transfTag;
  if (OPS_GetIntInput(&numData, &transfTag) < 0) {
    opserr << "WARNING transfTag is not integer\n";
    return 0;
  }

  // options
  double HirI = 0.0, HirJ = 0.0, HkrA = 0.0, HkrI = 0.0, HkrJ = 0.0;
  double mass = 0.0, Wtol = 1e-16, yftol = 1e-8;
  int cMass = 0, MaxIter = 15, nrow = 4;
  double coefData[3 * 4];
  for (int i = 0; i < 3 * nrow; i++) {
    coefData[i] = 0.0;
  }
  coefData[0] = 1.0;
  coefData[1] = 1.0;
  coefData[2] = 3.5;
  coefData[3] = -1.0;
  coefData[4] = 2.0;
  coefData[6] = 2.0;
  coefData[9] = 2.0;
  coefData[10] = 2.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if (type == "-Wtol") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &Wtol) < 0)
          return 0;
      }
    } else if (type == "-yftol") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &yftol) < 0)
          return 0;
      }
    } else if (type == "-MaxIter") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &MaxIter) < 0)
          return 0;
      }
    } else if (type == "-HirI") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &HirI) < 0)
          return 0;
      }
    } else if (type == "-HirJ") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &HirJ) < 0)
          return 0;
      }
    } else if (type == "-HkrA") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &HkrA) < 0)
          return 0;
      }
    } else if (type == "-HkrI") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &HkrI) < 0)
          return 0;
      }
    } else if (type == "-HkrJ") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &HkrJ) < 0)
          return 0;
      }
    } else if (type == "-mass") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0)
          return 0;
      }
    } else if (type == "-cMass") {
      cMass = 1;
    } else if (type == "-ySurf") {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &nrow) < 0)
          return 0;
        numData = 3 * nrow;
        if (OPS_GetDoubleInput(&numData, &coefData[0]) < 0)
          return 0;
        numData = 1;
      }
    }
  }

  // Create vector for isotropic hardening ratio
  Vector Hir(2);
  Hir(0) = HirI;
  Hir(1) = HirJ;

  // Create vector for kinematic hardening ratio
  Vector Hkr(3);
  Hkr(0) = HkrA;
  Hkr(1) = HkrI;
  Hkr(2) = HkrJ;

  // Create matrix with coefficients of yield surface
  Matrix GPYSC(coefData, nrow, 3);

  // check transf
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "coord transformation not found\n";
    return 0;
  }

  return new BeamColumnwLHNMYS(iData[0], iData[1], iData[2], data[0], data[1],
                               data[2], data[3], data[4], data[5], data[6],
                               *theTransf, yftol, Wtol, MaxIter, Hir, Hkr, nrow,
                               GPYSC, mass, cMass);
}
