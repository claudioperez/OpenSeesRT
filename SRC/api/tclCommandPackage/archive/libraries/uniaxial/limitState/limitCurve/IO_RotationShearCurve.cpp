

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/limitCurve/RotationShearCurve.h>
void *OPS_RotationShearCurve(void)
{
  if (shearCurveCount == 0) {
    // opserr << "RotationShearCurve limit curve - Written by MRL UT Austin
    // Copyright 2012 -  Use at your Own Peril \n";
    shearCurveCount++;
  }

  int argc = OPS_GetNumRemainingInputArgs();

  if (!(argc == 9 || argc == 23)) {
    opserr << "WARNING RotationShearCurve -- insufficient arguments\n";
    opserr << "For direct input of shear curve parameters and degrading slope "
              "want:\n\n";
    opserr << "limitCurve RotationShearCurve crvTag? eleTag? \n";
    opserr << "ndI? ndJ? rotAxis? Vn? Vr? Kdeg? rotLim? \n" << endln;

    opserr << "OR for calibrated shear curve and degrading slope want:\n\n";
    opserr << "limitCurve RotationShearCurve crvTag? eleTag?\n";
    opserr << "ndI? ndJ? rotAxis? Vn? Vr? Kdeg? defType?\n";
    opserr << "b? d? h? L? st? As? Acc? ld? db? rhot? f'c?\n";
    opserr << "fy? fyt? delta?\n" << endln;

    return 0;
  }

  int iTagData[2];
  int iNodeData[3];
  double dKdegData[3];
  double dRotLimData[1];
  int iTypeData[1];
  double dPropData[14];

  int numData;
  numData = 2;
  if (OPS_GetIntInput(&numData, iTagData) != 0) {
    opserr << "WARNING RotationShearCurve -- invalid crvTag? eleTag?\n"
           << endln;
    return 0;
  }
  int eleTag = iTagData[1];
  Domain *theDomain = 0;
  theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr
        << "WARNING RotationShearCurve -- Pointer to Domain was not returned\n"
        << endln;
    return 0;
  }
  Element *theElement = 0;
  theElement = theDomain->getElement(eleTag);
  if (theElement == 0) {
    opserr << "WARNING RotationShearCurve -- Element with tag " << iTagData[1]
           << " does not exist for shear curve tag " << iTagData[0] << endln
           << endln;
    return 0;
  }
  numData = 3;
  if (OPS_GetIntInput(&numData, iNodeData) != 0) {
    opserr << "WARNING RotationShearCurve -- invalid ndI? ndJ? rotAxis?\n"
           << endln;
    return 0;
  }
  Node *theNodeI = 0;
  theNodeI = theDomain->getNode(iNodeData[0]);
  if (theNodeI == 0) {
    opserr << "WARNING RotationShearCurve -- Node with tag " << iNodeData[0]
           << " does not exist for shear curve tag " << iTagData[0] << endln
           << endln;
    return 0;
  }
  Node *theNodeJ = 0;
  theNodeJ = theDomain->getNode(iNodeData[1]);
  if (theNodeJ == 0) {
    opserr << "WARNING RotationShearCurve -- Node with tag " << iNodeData[1]
           << " does not exist for shear curve tag " << iTagData[0] << endln
           << endln;
    return 0;
  }
  if (iNodeData[2] < 3 || iNodeData[2] > 6) {
    opserr << "WARNING RotationShearCurve -- rotAxis is invalid\n";
    opserr << "rotAxis = 3 -- Rotation about z-axis - 2D\n";
    opserr << "rotAxis = 4 -- Rotation about x-axis - 3D\n";
    opserr << "rotAxis = 5 -- Rotation about y-axis - 3D\n";
    opserr << "rotAxis = 6 -- Rotation about z-axis - 3D\n" << endln;
    return 0;
  }
  if (argc == 9) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dKdegData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid Vn? Vr? Kdeg?\n"
             << endln;
      return 0;
    }
    if (dKdegData[0] != -1 && !(dKdegData[0] > 0)) {
      opserr << "WARNING RotationShearCurve --  Vn input is invalid\n";
      opserr << "Vn = -1 -- Shear critical limit is not used\n";
      opserr << "Vn > 0 -- Shear critical limit is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[1] < -1) {
      opserr << "WARNING RotationShearCurve -- Vr input is invalid\n";
      opserr << "Vr = -1 -- Residual shear strength = 0.2*(maximum shear at "
                "failure)\n";
      opserr << "-1 < Vr < 0 -- Residual shear strength = Vr*(maximum shear at "
                "failure)\n";
      opserr << "Vr >= 0 -- Residual shear strength is the input value\n"
             << endln;
      return 0;
    }
    if (dKdegData[2] >= 0) {
      opserr << "WARNING RotationShearCurve -- Kdeg input is invalid\n";
      opserr << "The degrading slope must be less than zero\n" << endln;
      return 0;
    }
    numData = 1;
    if (OPS_GetDoubleInput(&numData, dRotLimData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid rotLim?\n" << endln;
      return 0;
    }
    if (dRotLimData[0] <= 0) {
      opserr << "WARNING RotationShearCurve -- rotLim input must be greater "
                "than zero\n"
             << endln;
      return 0;
    }

    LimitCurve *theCurve = 0;
    theCurve = new RotationShearCurve(
        iTagData[0], iTagData[1], iNodeData[0], iNodeData[1], iNodeData[2],
        dKdegData[0], dKdegData[1], dKdegData[2], dRotLimData[0], 0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, theDomain,
        theElement, theNodeI, theNodeJ);
    if (theCurve == 0) {
      opserr << "WARNING RotationShearCurve -- could not create limitCurve "
                "with constructor "
             << iTagData[0] << endln << endln;
      return 0;
    }
    return theCurve;
  } else {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dKdegData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid Vn? Vr? Kdeg?\n"
             << endln;
      return 0;
    }
    if (dKdegData[0] != -1 && !(dKdegData[0] >= 0)) {
      opserr << "WARNING RotationShearCurve --  Vn input is invalid\n";
      opserr << "Vn = -1 -- Shear critical limit is not used\n";
      opserr << "Vn = 0 -- Shear critical limit is calculated using ASCE 41 "
                "Eq. 6-4\n";
      opserr << "Vn > 0 -- Shear critical limit is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[1] < -1) {
      opserr << "WARNING RotationShearCurve -- Vr input is invalid\n";
      opserr << "Vr = -1 -- Residual shear strength from regression\n";
      opserr << "-1 < Vr < 0 -- Residual shear strength = Vr*(maximum shear at "
                "failure)\n";
      opserr << "Vr >= 0 -- Residual shear strength is the input value\n"
             << endln;
      return 0;
    }
    if (dKdegData[2] > 0) {
      opserr << "WARNING RotationShearCurve -- Kdeg input is invalid\n";
      opserr << "Kdeg = 0 -- Degrading slope calculated by regressions\n";
      opserr << "Kdeg < 0 -- Degrading slope is the input value\n" << endln;
      return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, iTypeData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid defType?\n" << endln;
      return 0;
    }
    if (iTypeData[0] > 5 || iTypeData[0] <= 0) {
      opserr << "WARNING RotationShearCurve -- invalid defType input?\n"
             << endln;
      opserr
          << "1 -- Flexure-Shear capacity based on theta_f rotation capacity\n";
      opserr << "2 -- Flexure-Shear capacity based on theta_total rotation "
                "capacity\n";
      opserr << "3 -- Flexure-Shear capacity based on theta_flexural rotation "
                "capacity\n";
      opserr << "4 -- Flexure-Shear capacity based on theta_total-plastic "
                "rotation capacity\n";
      opserr << "5 -- Flexure-Shear capacity based on theta_flexural-plastic "
                "rotation capacity\n"
             << endln;
      return 0;
    }
    numData = 14;
    if (OPS_GetDoubleInput(&numData, dPropData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid b? d? h? L? st? As? "
                "Acc? ld? db? rhot? f'c? fy? fyt? delta?\n"
             << endln;
      return 0;
    }
    LimitCurve *theCurve = 0;
    theCurve = new RotationShearCurve(
        iTagData[0], iTagData[1], iNodeData[0], iNodeData[1], iNodeData[2],
        dKdegData[0], dKdegData[1], dKdegData[2], 0.0, iTypeData[0],
        fabs(dPropData[0]), fabs(dPropData[1]), fabs(dPropData[2]),
        fabs(dPropData[3]), fabs(dPropData[4]), fabs(dPropData[5]),
        fabs(dPropData[6]), fabs(dPropData[7]), fabs(dPropData[8]),
        fabs(dPropData[9]), fabs(dPropData[10]), fabs(dPropData[11]),
        fabs(dPropData[12]), dPropData[13], theDomain, theElement, theNodeI,
        theNodeJ);
    if (theCurve == 0) {
      opserr << "WARNING RotationShearCurve -- could not create limitCurve "
                "with constructor "
             << iTagData[0] << endln << endln;
      return 0;
    }
    return theCurve;
  }
}
