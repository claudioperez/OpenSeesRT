

#include <g3_api.h>


#include <SRC/material/uniaxial/ConfinedConcrete01.h>
OPS_Export void *OPS_ConfinedConcrete01Material()
{
  if (numConfinedConcrete01Materials == 0) {
    numConfinedConcrete01Materials++;
    opserr << "ConfinedConceret01 unaxial material - Written by M.D'Amato, "
              "University of Basilicata, Italy 2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  /*
    USER INPUT PARAMETERS

    ----Unconfined Concrete Properties----
    1.  fc
    3.  epscuOption:
       0-direct input of epscu
       1-Scott
       2-gamma
    4.  epscu: input dependent upon epscuOption
    5.  nuOption: 0-constant
                          1-variable, upper bound 0.5
                          2-variable, no upper bound
    6.  nuc: poisson's ratio (if nuOption is 0)

    ---Section Geometric Properties---
    7.  secType: see Fig. 8 in the paper (Journal of Structural Engineering,
    page 1408)
    8.  dim: the number of square sections (determined from secType)


    ----Longitudinal Reinforcement Properties----
    17. phiLon: longitudinal bar diameter

    ----Attard and Setunge Properties----
    18. aggrType: 0-crushed
                  1-gravel
    19. concrType: 0-without silica fume
                   1-with silica fume

    ----Transverse Reinforcement Properties----

    9.  semiLength[]: semilengths of the square section
    10. phis[]: transverse bar diameter or wrapping thickness
    11. S[] : spacing of transverse reinforcement
    12. fyh[]: yield strength of transverse reinforcement
    13. mueps[]: ductility factor of transverse reinforcement
    14. Es0[]: initial elastic modulus of transverse reinforcement
    15. haRatio[]: hardening ratio of transverse reinforcement

    ---TCL Command Reference----------------------------

    OLD:
    uniaxialMaterial ConfinedConcrete01 $tag $secType $fpc <-stRatio $stRatio>
    $Ec
    <$epscu> <-scott> <-gamma $gamma <-epscu $epscuLimit>>   ;one option
    required
    <$nu> <-varub> <-varnoub>	;one option required
    $L1 <$L2> <$L3>	  ;Geometric properties
    $phis $S $fyh $Es0 $haRatio <-mu $mu>   ;External transverse
    <-phi $phisi> <-S $Si> <-fyh $fyhi> <-Es0 $Es0i> <-haRatio $haRatioi> <-mu
    $mui>	;Internal transverse
    <-wrap $cover $Am $Sw $ful $Es0w>  ;External wrapping
    $phiLon
    <-gravel> <-silica> <-tol $tol> <-maxNumIter $maxNumIter> ;Optional

    PROPOSED:
    uniaxialMaterial ConfinedConcrete01 $tag $secType $fpc $Ec
    <-epscu $epscu> <-scott> <-gamma $gamma>  ;one option required
    <-nu $nu> <-varub> <-varnoub> ;one option required
    $L1 ($L2) ($L3) ; provision of L2 and L3 depend on previous args passed
    $phis $S $fyh $Es0 $haRatio $mu
    $phiLon
    FOLLOWING ARE OPTIONAL
    <-internal $phisi $Si $fyhi $Es0i $haRatioi $mui
    <-wrap $cover $Am $Sw $ful $Es0w>
    <-gravel> <-silica> <-tol $tol> <-maxNumIter $maxNumIter>  <-epscuLimit
    $epscuLimit>
    <-stRatio $stRatio>

*/
  int argLoc = 4;
  int numReq = 12;
  double temp;
  int tag;
  int secType;
  int dim = 1;
  std::vector<double> semiLength, phis, S, fyh, Es0, haRatio, mueps, As, Is;
  double rhos;
  double fpc;
  double stRatio = 1.0;
  double Ec;
  int epscuOption; // 0-direct input of epscu
  // 1-Scott
  // 2-gamma
  double epscu;     // input dependent upon epscuOption
  int nuOption = 0; // 0-constant
  // 1-variable, upper bound 0.5
  // 2-variable, no upper bound
  double nuc = 0; // poisson's ratio (if nuOption is 0)
  double phiLon;
  int concrType = 0;
  int aggrType = 0;

  double tol = 0.000001;
  int maxNumIter = 500;

  int argc = OPS_GetNumRemainingInputArgs();

  if (argc < numReq) {
    opserr << "WARNING insufficient arguments\n";
    return 0;
  }

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag"
           << endln;
    return 0;
  }

  // ==Read required ConfinedConcrete01 material
  // parameters===========================

  // --Parse section type-----------------------------------------------------
  //  char argvS[5];
  const char *argvS;
  argvS = OPS_GetString();
  /*
  if (OPS_GetString(argvS, 5) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" <<
  endln; return 0;
  }
  */

  if (strcmp(argvS, "S1") == 0 || strcmp(argvS, "s1") == 0) {
    secType = 1;
  } else if (strcmp(argvS, "S2") == 0 || strcmp(argvS, "s2") == 0) {
    secType = 2;
  } else if (strcmp(argvS, "S3") == 0 || strcmp(argvS, "s3") == 0) {
    secType = 3;
  } else if (strcmp(argvS, "S4a") == 0 || strcmp(argvS, "s4a") == 0) {
    secType = 41;
  } else if (strcmp(argvS, "S4b") == 0 || strcmp(argvS, "s4b") == 0) {
    secType = 42;
  } else if (strcmp(argvS, "S5") == 0 || strcmp(argvS, "s5") == 0) {
    secType = 5;
  } else if (strcmp(argvS, "C") == 0 || strcmp(argvS, "c") == 0) {
    secType = 6;
  } else if (strcmp(argvS, "R") == 0 || strcmp(argvS, "r") == 0) {
    secType = 7;
  } else {
    opserr << "WARNING invalid section type, should be:(S1-S3, S4a, S4b, S5, "
              "R, C)\n";
    opserr << "ConfinedConcrete01 material: " << tag << endln;
    return 0;
  }

  // --Parse concrete
  // properties-------------------------------------------------
  if (OPS_GetDouble(&numData, &fpc) != 0) {
    opserr << "WARNING invalid fpc\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }

  if (fpc < 0.0)
    fpc = -fpc;

  // printf("fpc: %f\n", fpc);

  if (OPS_GetDouble(&numData, &Ec) != 0) {
    opserr << "WARNING invalid Ec\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }

  // printf("Ec: %f\n", Ec);

  // --Parse epsilon_cu: 3
  // options-----------------------------------------------
  double epscuLimit = 0.05;

  // char argvEPSCU[10];
  const char *argvEPSCU = OPS_GetString();
  argLoc++;

  /*
  if (OPS_GetString(argvEPSCU, 10) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" <<
  endln; return 0; } else

  */

  if (strcmp(argvEPSCU, "-scott") == 0 || strcmp(argvEPSCU, "-Scott") == 0) {
    epscuOption = 1;
  } else if (strcmp(argvEPSCU, "-gamma") == 0) {
    // Update number of arguments required
    numReq++;
    if (argc < numReq) {
      opserr << "WARNING insufficient arguments\n";
      return 0;
    }
    epscuOption = 2;
    if (OPS_GetDouble(&numData, &epscu) != 0) {
      opserr << "WARNING invalid gamma\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    } else
      argLoc++;

  } else { // Specify epscu directly
    epscuOption = 0;
    if (OPS_GetDouble(&numData, &epscu) != 0) {
      opserr << "WARNING invalid epscu\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    } else
      argLoc++;
  }

  if (epscuOption != 1 && epscu < 0.0)
    epscu = -epscu;

  // --Parse nu----------------------------------------------------------------

  //  char argvNu[10];
  const char *argvNu = OPS_GetString();
  /*
  if (OPS_GetString(argvNu, 10) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" <<
  endln; return 0; } else

  */
  argLoc++;

  if (strcmp(argvNu, "-varUB") == 0 || strcmp(argvNu, "-varub") == 0) {
    nuOption = 1;
  } else if (strcmp(argvNu, "-varNoUB") == 0 ||
             strcmp(argvNu, "-varnoub") == 0) {
    nuOption = 2;
  } else {
    nuOption = 0;
    if (OPS_GetDouble(&numData, &nuc) != 0) {
      opserr << "WARNING invalid nu\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    } else
      argLoc++;
  }

  // --Parse section geometry--------------------------------------------------
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid L1\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  } else
    argLoc++;

  semiLength.push_back(temp / 2);

  if (secType == 2 || secType == 3) {
    dim++;
    semiLength.push_back(semiLength[0] * sqrt(2.0) / 2);
  } else if (secType == 5) {
    dim++;
    semiLength.push_back(semiLength[0]);
  }

  if (secType == 42) {
    dim++;
    semiLength.push_back(semiLength[0]);
  }

  if (secType == 41 || secType == 42 || secType == 7) {
    dim++;
    if (OPS_GetDouble(&numData, &temp) != 0) {
      opserr << "WARNING invalid L2\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    }
    semiLength.push_back(temp / 2);
  }

  if (secType == 41) {
    dim += 2;
    semiLength.push_back(semiLength[0]); // semiLength[2] = semiLength[0]
    if (OPS_GetDouble(&numData, &temp) != 0) {
      opserr << "WARNING invalid L3\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    }
    semiLength.push_back(temp / 2);
  }

  // --Parse transverse reinforcement
  // properties---------------------------------
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid phis\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  phis.push_back(temp);
  As.push_back(PI * pow(phis[0], 2) / 4);
  Is.push_back(PI * pow(phis[0], 4) / 64);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid S\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  S.push_back(temp);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid fyh\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  fyh.push_back(temp);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid Es0\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  Es0.push_back(temp);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid haRatio\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  haRatio.push_back(temp);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid mu\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  mueps.push_back(temp);

  if (OPS_GetDouble(&numData, &phiLon) != 0) {
    opserr << "WARNING invalid phiLon \n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }

  // Preload transverse properties in other dimensions
  for (int i = 1; i < dim; i++) {
    if (secType == 2) {
      As.push_back(As[0] * sqrt(2.0) / 2.0 / 2.0);
      Is.push_back(0.0);
    } else {
      As.push_back(As[0]);
      Is.push_back(Is[0]);
    }
    phis.push_back(phis[0]);
    S.push_back(S[0]);
    fyh.push_back(fyh[0]);
    mueps.push_back(mueps[0]);
    Es0.push_back(Es0[0]);
    haRatio.push_back(haRatio[0]);
  }

  // LATER
  argc = OPS_GetNumRemainingInputArgs();
  while (argc > 0) {
    const char *argvLoc = OPS_GetString();
    ;
    /*
    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" <<
    endln; return 0;
    }
    */

    if (strcmp(argvLoc, "-stRatio") == 0) {
      if (OPS_GetDouble(&numData, &stRatio) != 0 || stRatio > 1.0 ||
          stRatio < 0.0) {
        opserr << "WARNING invalid stRatio\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
    }

    else if (strcmp(argvLoc, "-epscuLimit") == 0) {
      if (OPS_GetDouble(&numData, &epscuLimit) != 0) {
        opserr << "WARNING invalid epscuLimit\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
    }

    // internal transverse
    else if ((strcmp(argvLoc, "-internalT") == 0) ||
             (strcmp(argvLoc, "-internal") == 0)) {
      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid phi (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      phis[dim - 1] = temp;
      if (secType == 2) {
        As[dim - 1] = temp * sqrt(2.0) / 2.0 / 2.0;
        Is[dim - 1] = 0.0;
      } else {
        As[dim - 1] = PI * pow(temp, 2) / 4;
        Is[dim - 1] = PI * pow(temp, 4) / 64;
      }

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid S (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      S[dim - 1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid fyh (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      fyh[dim - 1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid Es0 (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      Es0[dim - 1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid haRatio (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      haRatio[dim - 1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid mu (stirrups)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      mueps[dim - 1] = temp;

      // One internal transverse property parsed
      // Check wrapping for allowable sections
      if (!((2 <= secType && secType <= 5) || secType == 41 || secType == 42)) {
        opserr << "WARNING Stirrups only valid for section types S2, S3, S4a, "
                  "S4b, S5\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }

      if (secType == 7) {
        phis[dim - 2] = phis[dim - 1];
        As[dim - 2] = As[dim - 1];
        Is[dim - 2] = Is[dim - 1];
        S[dim - 2] = S[dim - 1];
        fyh[dim - 2] = fyh[dim - 1];
        Es0[dim - 2] = Es0[dim - 1];
        haRatio[dim - 2] = haRatio[dim - 1];
        mueps[dim - 2] = mueps[dim - 1];
      }
    }

    else if (strcmp(argvLoc, "-wrap") == 0) {
      if (argc < 6) {
        opserr << "WARNING insufficient arguments with -wrap option\n";
        return 0;
      }

      dim++;
      phis.push_back(0.0);
      Is.push_back(0.0);
      mueps.push_back(1.0);
      haRatio.push_back(1.0);
      argLoc++;

      // Get the cover, update semiLength vector
      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid cover\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      semiLength.push_back(semiLength[0] + phis[0] / 2 + temp);

      if (secType == 7) {
        dim++;
        semiLength.push_back(semiLength[1] + phis[1] / 2 + temp);
      }

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid As (wrapping)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      As.push_back(temp);

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid S (wrapping) \n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      S.push_back(temp);

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid fyh (wrapping)\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      fyh.push_back(temp);

      if (OPS_GetDouble(&numData, &temp) != 0) {
        opserr << "WARNING invalid Es0 (wrapping) \n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }
      Es0.push_back(temp);

      if (secType == 7) {
        phis.push_back(phis[dim - 2]);
        As.push_back(As[dim - 2]);
        Is.push_back(Is[dim - 2]);
        S.push_back(S[dim - 2]);
        fyh.push_back(fyh[dim - 2]);
        mueps.push_back(mueps[dim - 2]);
        Es0.push_back(Es0[dim - 2]);
        haRatio.push_back(haRatio[dim - 2]);
      }
    }

    else if (strcmp(argvLoc, "-gravel") == 0) {
      aggrType = 1;

    } else if (strcmp(argvLoc, "-silica") == 0) {
      concrType = 1;

    } else if (strcmp(argvLoc, "-tol") == 0) {

      if (OPS_GetDouble(&numData, &tol) != 0) {
        opserr << "WARNING invalid tol\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }

    } else if (strcmp(argvLoc, "-maxNumIter") == 0) {
      if (OPS_GetInt(&numData, &maxNumIter) != 0) {
        opserr << "WARNING invalid maxNumIter\n";
        opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
        return 0;
      }

    } else {
      opserr << "WARNING invalid argument(s) :" << argvLoc << endln;
      return 0;
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  //--Transverse ratio--------------------------------------------------------
  if (secType == 1) { // S1
    rhos = 2 * As[0] /
           (S[0] * semiLength[0]); // Add. External wrapping is neglected
  } else if (secType == 2) {       // S2
    rhos = (2 * As[0] + As[1]) / (S[0] * semiLength[0]);
  } else if (secType == 3) { // S3
    rhos = (2 * As[0] + sqrt(2.0) * As[1]) / (S[0] * semiLength[0]);
  } else if (secType == 41) { // S41
    rhos = (3 * As[0] + As[1]) / (S[0] * semiLength[0]);
  } else if (secType == 42) { // S42
    rhos = 2 * (As[0] + As[1]) / (S[0] * semiLength[0]);
  } else if (secType == 5) { // S5
    rhos = 2 * (As[0] + (1 + sqrt(2.0)) / 3.0 * As[1]) / (S[0] * semiLength[0]);
  } else if (secType == 6) { // S6
    rhos = 2 * As[0] / (S[0] * semiLength[0]);
    // Add. External wrapping is neglected
  } else { // S7
    rhos = As[0] * (semiLength[0] + semiLength[1]) /
           (S[0] * semiLength[0] * semiLength[1]);
    // Add. External wrapping is neglected
  }

  // printf("rhos %f\n", rhos);

  // Parse additional optional parameters

  if (epscuOption == 1) {
    epscu = 0.004 + 0.9 * rhos * (fyh[0] / 300.0); // Scott et al. 1982
  }

  /*
    printf("SemiLengths: \n");

    for (int n = 0; n < semiLength.size(); n++) {
    printf("%f\n", semiLength[n]);
    }

    printf("Dimension: %d\n", dim);
    printf("Nu: %f\n", nuc);
  */

  // Parsing was successful, allocate the material
  theMaterial = new ConfinedConcrete01(
      tag, secType, dim, semiLength, phis, S, fyh, Es0, haRatio, mueps, As, Is,
      rhos, fpc, stRatio, Ec, epscuOption, epscu, epscuLimit, nuOption, nuc,
      phiLon, concrType, aggrType, tol, maxNumIter);

  return theMaterial;
}
