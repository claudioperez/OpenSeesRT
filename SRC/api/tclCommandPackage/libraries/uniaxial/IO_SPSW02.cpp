

#include <g3_api.h>


#include <SRC/material/uniaxial/SPSW02.h>
void *OPS_SPSW02(G3_Runtime *rt)
{
  if (numSPSWcall == 0) {
    opserr << "------ SPSW02 unaxialMaterial, Written by SAJalali @ Amirkabir "
              "University of Technology, Tehran, 2015-------\n";
    opserr << "------------------------------ Please Send Comments to: "
              "seyyed-jalali@aut.ac.ir-----------------------------\n";
    opserr << "-------Syntax:\n";
    opserr << "-------UniaxialMaterial SPSW02 tag ";
    opserr << "-------E0 b <-geom Fpy t h l> <-params Fts Fcs cmpUnldngEFac "
              "sigTEFac sigTFfac epsTFfac> -R $R -Damage epsPCFac pstCapEFac "
              "gama c resFac\n\n";
    opserr << "----------------------------------------------------------------"
              "--------------------------------------------\n\n\n";
    numSPSWcall = 1;
  }
  int tag;
  double fpy, E0, b, t, hs, R, l, gama, c, epsPCFac, pstCapEFac, resFac;
  double Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac;
  fpy = 0;
  bool paramsSet = false;
  UniaxialMaterial *theMaterial = 0;
  int argc = OPS_GetNumRemainingInputArgs();
  int numData = 1;
  int curArg = 2;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid -tag";
    opserr << "uniaxialMaterial SPSW02: " << tag << endln;
    return 0;
  }
  curArg++;

  if (OPS_GetDoubleInput(&numData, &E0) != 0) {
    opserr << "WARNING invalid -E0";
    opserr << "uniaxialMaterial SPSW02: " << tag << endln;
    return 0;
  }
  curArg++;

  if (OPS_GetDoubleInput(&numData, &b) != 0) {
    opserr << "WARNING invalid -b";
    opserr << "uniaxialMaterial SPSW02: " << tag << endln;
    return 0;
  }
  curArg++;

  const char *str = OPS_GetString();
  curArg++;

  if (strcmp(str, "-geom") == 0) {
    if (OPS_GetDoubleInput(&numData, &fpy) != 0) {
      opserr << "WARNING invalid -Fts";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &t) != 0) {
      opserr << "WARNING invalid -t";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &hs) != 0) {
      opserr << "WARNING invalid -h";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &l) != 0) {
      opserr << "WARNING invalid -l";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;
  } else if (strcmp(str, "-params") == 0) {
    paramsSet = true;
    // Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac, epsTFfac
    if (OPS_GetDoubleInput(&numData, &Fts) != 0) {
      opserr << "WARNING invalid Fts";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &Fcs) != 0) {
      opserr << "WARNING invalid Fcs";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &cmpUnldngEFac) != 0) {
      opserr << "WARNING invalid cmpUnldngEFac";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &sigTEFac) != 0) {
      opserr << "WARNING invalid sigTEFac";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &sigTFfac) != 0) {
      opserr << "WARNING invalid sigTFfac";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;

    if (OPS_GetDoubleInput(&numData, &epsTFfac) != 0) {
      opserr << "WARNING invalid epsTFfac";
      opserr << "uniaxialMaterial SPSW02: " << tag << endln;
      return 0;
    }
    curArg++;
  }
  if (fpy == 0 && !paramsSet) {
    opserr
        << "WARNING at least one of -params or -geom options must be provided";
    opserr << "uniaxialMaterial SPSW02: " << tag << endln;
    return 0;
  }
  if (fpy != 0 && paramsSet) {
    opserr << "WARNING both -params and -geom options cannot be used at the "
              "same time";
    opserr << "uniaxialMaterial SPSW02: " << tag << endln;
    return 0;
  }
  R = 50;

  if (argc >= curArg + 1) {
    /*char[] str[10];
    OPS_GetString(str, 5);*/		//for OpenSees 2.4.5 and older
    str = OPS_GetString();
    curArg++;
    if (strcmp(str, "-R") == 0) {
      if (OPS_GetDoubleInput(&numData, &R) != 0) {
        opserr << "WARNING invalid -R";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;
    }
  }
  epsPCFac = 1.e20;
  pstCapEFac = b;
  gama = 10000;
  c = 1;
  resFac = 0.001;
  if (argc >= curArg + 1) {
    /*char[] str[10];
    OPS_GetString(str, 124);*/
    str = OPS_GetString();
    curArg++;
    if (strcmp(str, "-Damage") == 0 || strcmp(str, "-damage") == 0) {
      if (OPS_GetDoubleInput(&numData, &epsPCFac) != 0) {
        opserr << "WARNING invalid -epsPCFac";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;

      if (OPS_GetDoubleInput(&numData, &pstCapEFac) != 0) {
        opserr << "WARNING invalid -pstCapEFac";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;

      if (OPS_GetDoubleInput(&numData, &gama) != 0) {
        opserr << "WARNING invalid -gama";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;

      if (OPS_GetDoubleInput(&numData, &c) != 0) {
        opserr << "WARNING invalid -c";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;

      if (OPS_GetDoubleInput(&numData, &resFac) != 0) {
        opserr << "WARNING invalid -resFac";
        opserr << "uniaxialMaterial SPSW02: " << tag << endln;
        return 0;
      }
      curArg++;
    }
  }
  if (paramsSet)
    theMaterial =
        new SPSW02(tag, E0, b, Fts, Fcs, cmpUnldngEFac, sigTEFac, sigTFfac,
                   epsTFfac, R, epsPCFac, pstCapEFac, gama, c, resFac);
  else
    theMaterial = new SPSW02(tag, fpy, E0, b, t, hs, l, R, epsPCFac, pstCapEFac,
                             gama, c, resFac);
  // opserr<<"ok\n";

  return theMaterial;
}
