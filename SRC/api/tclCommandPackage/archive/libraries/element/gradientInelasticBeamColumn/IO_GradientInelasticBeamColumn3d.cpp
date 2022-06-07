

#include <g3_api.h>


#include <SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.h>
void *OPS_GradientInelasticBeamColumn3d()
{
  // Necessary Arguments
  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr
        << "WARNING! gradientInelasticBeamColumn3d - insufficient arguments\n"
        << "         Want: eleTag? iNode? jNode? transfTag? integrationTag? "
           "lambda1? lambda2? lc?\n"
        << "         <-constH> <-iter maxIter? minTol? maxTol?> <-corControl "
           "maxEpsInc? maxPhiInc?>\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - ndm must be 3 and ndf "
              "must be 6\n";
    return 0;
  }

  // inputs:
  int iData[5];
  int numData = 5;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - invalid input tags\n";
    return 0;
  }

  int eleTag = iData[0];
  int nodeTagI = iData[1];
  int nodeTagJ = iData[2];
  int transfTag = iData[3];
  int integrTag = iData[4];

  double ddata[3];
  numData = 3;
  if (OPS_GetDoubleInput(&numData, ddata) < 0) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - invalid lc\n";
    return 0;
  }
  double lam1 = ddata[0];
  double lam2 = ddata[1];
  double lc = ddata[2];

  // Optional Arguments
  int maxIter = 50;
  double minTol = 1E-10, maxTol = 1E-8;
  bool correctionControl = false;
  bool constH = false;
  double maxEpsInc = 0.0, maxPhiInc = 0.0;

  numData = 1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *word = OPS_GetString();

    if (strcmp(word, "-constH") == 0)
      constH = true;
    else if (strcmp(word, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 2) {
        if (OPS_GetIntInput(&numData, &maxIter) < 0) {
          opserr
              << "WARNING! gradientInelasticBeamColumn3d - invalid maxIter\n";
          return 0;
        }
        if (OPS_GetDoubleInput(&numData, &minTol) < 0) {
          opserr << "WARNING! gradientInelasticBeamColumn3d - invalid minTol\n";
          return 0;
        }
        if (OPS_GetDoubleInput(&numData, &maxTol) < 0) {
          opserr << "WARNING! gradientInelasticBeamColumn3d - invalid maxTol\n";
          return 0;
        }
      } else {
        opserr << "WARNING! gradientInelasticBeamColumn3d - need maxIter? "
                  "minTol? maxTol? after -iter \n";
        return 0;
      }
    } else if (strcmp(word, "-corControl") == 0) {
      correctionControl = true;

      if (OPS_GetNumRemainingInputArgs() > 1) {
        if (OPS_GetDoubleInput(&numData, &maxEpsInc) < 0) {
          opserr
              << "WARNING! gradientInelasticBeamColumn3d - invalid maxEpsInc\n";
          return 0;
        }
        if (OPS_GetDoubleInput(&numData, &maxPhiInc) < 0) {
          opserr
              << "WARNING! gradientInelasticBeamColumn3d - invalid maxPhiInc\n";
          return 0;
        }
      } else
        opserr << "WARNING! gradientInelasticBeamColumn3d - no max. correction "
                  "increments set\n"
               << "         -> setting them automatically|\n";
    }
  }

  // check transf
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - CrdTransf with tag "
           << transfTag << " not found\n";
    return 0;
  }

  // check beam integrataion
  BeamIntegrationRule *theRule = OPS_getBeamIntegrationRule(integrTag);
  if (theRule == 0) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - BeamIntegrationRule "
              "with tag "
           << integrTag << " not found\n";
    return 0;
  }

  BeamIntegration *beamIntegr = theRule->getBeamIntegration();
  if (beamIntegr == 0) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - failed to create beam "
              "integration\n";
    return 0;
  }

  // check sections
  const ID &secTags = theRule->getSectionTags();
  int numIntegrPoints = secTags.Size();

  for (int i = 2; i < numIntegrPoints; i++) {
    if (secTags(i) != secTags(i - 1)) {
      opserr << "WARNING! gradientInelasticBeamColumn3d - internal integration "
                "points should have identical tags\n"
             << "continued using section tag of integration point 2 for all "
                "internal integration points\n";
      return 0;
    }
  }

  SectionForceDeformation *endSection1 =
      OPS_getSectionForceDeformation(secTags(0));
  if (!endSection1) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - section with tag "
           << secTags(0) << " not found\n";
    return 0;
  }

  SectionForceDeformation *intSection =
      OPS_getSectionForceDeformation(secTags(1));
  if (!intSection) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - section with tag "
           << secTags(1) << " not found\n";
    return 0;
  }

  SectionForceDeformation *endSection2 =
      OPS_getSectionForceDeformation(secTags(numIntegrPoints - 1));
  if (!endSection2) {
    opserr << "WARNING! gradientInelasticBeamColumn3d - section with tag "
           << secTags(numIntegrPoints - 1) << " not found\n";
    return 0;
  }

  Element *theEle = new GradientInelasticBeamColumn3d(
      eleTag, nodeTagI, nodeTagJ, numIntegrPoints, &endSection1, &intSection,
      &endSection2, lam1, lam2, *beamIntegr, *theTransf, lc, minTol, maxTol,
      maxIter, constH, correctionControl, maxEpsInc, maxPhiInc);

  return theEle;
}
