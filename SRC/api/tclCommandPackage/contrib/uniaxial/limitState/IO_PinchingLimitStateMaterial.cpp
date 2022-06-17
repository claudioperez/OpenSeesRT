

#include <g3_api.h>


#include <SRC/material/uniaxial/limitState/PinchingLimitStateMaterial.h>
void *OPS_PinchingLimitState(void)
{
  if (numPinchingLimitStateMaterial == 0) {
    numPinchingLimitStateMaterial++;
    // opserr << "PinchingLimitStateMaterial unaxial material - Written by MRL
    // UT Austin Copyright 2012 - Use at Your Peril\n";
  }

  UniaxialMaterial *theMaterial = 0;
  LimitCurve *theCurve = 0;
  int argc = OPS_GetNumRemainingInputArgs();

  if (!(argc == 32 || argc == 21)) {
    opserr << "WARNING PinchingLimitStateMaterial -- insufficient arguments\n";
    opserr << "For direct input of limit state material want:\n\n";
    opserr << "uniaxialMaterial PinchingLimitStateMaterial matTag?\n";
    opserr << "nodeT? nodeB? driftAxis? Kelas? crvTyp? crvTag?\n";
    opserr << "YpinchUPN? YpinchRPN? XpinchRPN?\n";
    opserr << "YpinchUNP? YpinchRNP? XpinchRNP?\n";
    opserr << "dmgStrsLimE? dmgDispMax?\n?";
    opserr << "dmgE1? dmgE2? dmgE3? dmgE4? dmgELim?\n";
    opserr << "dmgR1? dmgR2? dmgR3? dmgR4? dmgRLim? dmgRCyc?\n";
    opserr << "dmgS1? dmgS2? dmgS3? dmgS4? dmgSLim? dmgSCyc?\n" << endln;

    opserr << "OR for calibrated limit state material want:\n\n";
    opserr << "uniaxialMaterial PinchingLimitStateMaterial matTag?\n";
    opserr << "nodeT? nodeB? driftAxis? Kelas? crvTyp? crvTag? eleTag?\n";
    opserr << "b? d? h? a? st? As? Acc? ld? db? rhot? f'c?\n";
    opserr << "fy? fyt?\n" << endln;

    return 0;
  }

  int iTagData[1];
  int iNodeData[3];
  double dKelasData[1];
  int iCrvData[2];
  double dpinchPN[3];
  double dpinchNP[3];
  double dDmgProp[2];
  double dDmgEdata[5];
  double dDmgRdata[6];
  double dDmgSdata[6];
  int iEleTag[1];
  double dPropData[13];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iTagData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid uniaxialMaterial "
              "matTag?\n"
           << endln;
    return 0;
  }
  numData = 3;
  if (OPS_GetIntInput(&numData, iNodeData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid nodeT? nodeB? "
              "driftAxis?\n"
           << endln;
    return 0;
  }
  Domain *theDomain = 0;
  theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- Pointer to Domain was not "
              "returned\n"
           << endln;
    return 0;
  }
  Node *theNodeT = 0;
  theNodeT = theDomain->getNode(iNodeData[0]);
  if (theNodeT == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- nodeT with tag "
           << iNodeData[0] << " does not exist for uniaxialMaterial tag "
           << iTagData[0] << endln << endln;
    return 0;
  }
  Node *theNodeB = 0;
  theNodeB = theDomain->getNode(iNodeData[1]);
  if (theNodeB == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- nodeB with tag "
           << iNodeData[1] << " does not exist for uniaxialMaterial tag "
           << iTagData[0] << endln << endln;
    return 0;
  }
  if (iNodeData[2] < 1 || iNodeData[2] > 3) {
    opserr << "WARNING PinchingLimitStateMaterial -- driftAxis is invalid\n";
    opserr << "driftAxis = 1 -- Drift along the x-axis\n";
    opserr << "driftAxis = 2 -- Drift along the y-axis\n";
    opserr << "driftAxis = 3 -- Drift along the z-axis\n";
    return 0;
  }
  numData = 1;
  if (OPS_GetDoubleInput(&numData, dKelasData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid Kelas?\n";
    return 0;
  }
  if ((dKelasData[0] < -4 || dKelasData[0] == 0) && argc == 23) {
    opserr << "WARNING PinchingLimitStateMaterial -- Kelas? is invalid\n";
    opserr << "Kelas = -4 -- Shear stiffness calculated assuming double "
              "curvature and shear springs top and bottom\n";
    opserr << "Kelas = -3 -- Shear stiffness calculated assuming double "
              "curvature and a shear spring at the bottom\n";
    opserr << "Kelas = -2 -- Shear stiffness calculated assuming single "
              "curvature and shear springs top and bottom\n";
    opserr << "Kelas = -1 -- Shear stiffness calculated assuming single "
              "curvature and a shear spring at the bottom\n";
    opserr << "Kelas > 0 -- Shear stiffness is the input value\n";
    return 0;
  } else if (dKelasData[0] <= 0 && argc == 34 /*39*/) {
    opserr << "WARNING PinchingLimitStateMaterial -- Kelas? is invalid\n";
    opserr << "Kelas must be greater than zero\n";
    return 0;
  }
  numData = 2;
  if (OPS_GetIntInput(&numData, iCrvData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid crvTyp? crvTag?\n"
           << endln;
    return 0;
  }
  int crvTyp = iCrvData[0];
  int crvTag = iCrvData[1];
  if (crvTyp == 2) {
    theCurve = OPS_getLimitCurve(crvTag);
    if (theCurve == 0 && crvTyp != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- limit curve with tag "
             << crvTag << " not found for material tag " << iTagData[0] << endln
             << endln;
      return 0;
    }
  }
  if (crvTyp < 0 || crvTyp > 2) {
    opserr << "WARNING PinchingLimitStateMaterial --  crvTyp? is invalid\n";
    opserr << "crvType = 0 -- no limit curve\n";
    opserr << "crvType = 1 -- axial limit curve\n";
    opserr << "crvType = 2 -- shear limit curve\n" << endln;
    return 0;
  }
  if (crvTyp == 1) {
    opserr << "WARNING PinchingLimitStateMaterial -- Axial curve has not been "
              "implemented\n"
           << endln;
    return 0;
  }
  if (argc == 32) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dpinchPN) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid YpinchUPN? "
                "YpinchRPN? XpinchRPN?\n"
             << endln;
      return 0;
    }
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dpinchNP) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid YpinchUNP? "
                "YpinchRNP? XpinchRNP?\n"
             << endln;
      return 0;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dDmgProp) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgStrsLimE? "
                "dmgDispMax?\n"
             << endln;
      return 0;
    }
    if (dDmgProp[0] < 0.0001)
      dDmgProp[0] = 0.0001;
    numData = 5;
    if (OPS_GetDoubleInput(&numData, dDmgEdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgE1? dmgE2? "
                "dmgE3? dmgE4? dmgELim?\n"
             << endln;
      return 0;
    }
    numData = 6;
    if (OPS_GetDoubleInput(&numData, dDmgRdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgR1? dmgR2? "
                "dmgR3? dmgR4? dmgRLim? dmgRCyc?\n"
             << endln;
      return 0;
    }
    numData = 6;
    if (OPS_GetDoubleInput(&numData, dDmgSdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgS1? dmgS2? "
                "dmgS3? dmgS4? dmgSLim? dmgSCyc?\n"
             << endln;
      return 0;
    }
    theMaterial = new PinchingLimitStateMaterial(
        iTagData[0], iNodeData[0], iNodeData[1], iNodeData[2], dKelasData[0],
        iCrvData[0], iCrvData[1], dpinchPN[0], dpinchPN[1], dpinchPN[2],
        dpinchNP[0], dpinchNP[1], dpinchNP[2], dDmgProp[0], dDmgProp[1],
        dDmgEdata[0], dDmgEdata[1], dDmgEdata[2], dDmgEdata[3], dDmgEdata[4],
        0.0, 0.0, 0.0, 0.0, 0.0, dDmgRdata[0], dDmgRdata[1], dDmgRdata[2],
        dDmgRdata[3], dDmgRdata[4], dDmgRdata[5], dDmgSdata[0], dDmgSdata[1],
        dDmgSdata[2], dDmgSdata[3], dDmgSdata[4], dDmgSdata[5], 0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, theDomain,
        theNodeT, theNodeB, *theCurve, 0);

    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial with "
                "PinchinLimitState\n";
      return 0;
    }
    return theMaterial;

  } else {
    numData = 1;
    if (OPS_GetIntInput(&numData, iEleTag) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid eleTag?\n"
             << endln;
      return 0;
    }
    int eleTag = iEleTag[0];
    Element *theElement = 0;
    theElement = theDomain->getElement(eleTag);
    if (theElement == 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- Element with tag "
             << eleTag << " does not exist for uniaxialMaterial tag "
             << iTagData[0] << endln << endln;
      return 0;
    }
    numData = 13;
    if (OPS_GetDoubleInput(&numData, dPropData) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid b? d? h? a? st? "
                "As? Acc? ld? db? rhot? f'c? fy? fyt?\n"
             << endln;
      return 0;
    }
    theMaterial = new PinchingLimitStateMaterial(
        iTagData[0], iNodeData[0], iNodeData[1], iNodeData[2], dKelasData[0],
        iCrvData[0], iCrvData[1], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, iEleTag[0], fabs(dPropData[0]),
        fabs(dPropData[1]), fabs(dPropData[2]), fabs(dPropData[3]),
        fabs(dPropData[4]), fabs(dPropData[5]), fabs(dPropData[6]),
        fabs(dPropData[7]), fabs(dPropData[8]), fabs(dPropData[9]),
        fabs(dPropData[10]), fabs(dPropData[11]), fabs(dPropData[12]),
        theDomain, theNodeT, theNodeB, *theCurve, theElement);

    if (theMaterial == 0) {
      opserr
          << "WARNING could not create uniaxialMaterial PinchingLimitState\n ";
      return 0;
    }
    return theMaterial;
  }
}
