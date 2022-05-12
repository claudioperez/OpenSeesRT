

#include <g3_api.h>


#include <SRC/material/uniaxial/Steel4.h>
void *OPS_Steel4(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double basicData[2], kinData[8], isoData[9], ultData[4], initData[1];
  int memData[1];
  int argc = 1, numBasic = 2, numKin = 4, numIso = 5, numUlt = 2, numMem = 1,
      numInit = 1;

  if (OPS_GetIntInput(&argc, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel4 tag" << endln;
    return 0;
  }

  // Read the basic parameters
  argc = OPS_GetNumRemainingInputArgs();
  // if the number of arguments is less than the minimum, throw an error
  if (argc < numBasic) {
    opserr << "Invalid #args, want at least two args for Steel4 in the "
              "following format:\n"
           << "uniaxialMaterial Steel4" << iData[0] << " E0? fy?" << endln;
    return 0;
  }
  // if the first two parameters are not doubles, throw an error
  else if (OPS_GetDoubleInput(&numBasic, basicData) != 0) {
    opserr << "Invalid args; E0 and fy for Steel4 (tag: " << iData[0]
           << ") shall be given as floating point numbers" << endln;
    return 0;
  }
  // Parsing was successful, the basic data are stored in basicData
  // Load the default values for everything else
  else {
    // kinematic hardening - bilinear perfectly plastic behavior by default
    kinData[0] = 0.0;  // b_k
    kinData[1] = 50.0; // R_0
    kinData[2] = 0.1;  // r_1
    kinData[3] = 0.15; // r_2
    for (unsigned int i = 0; i < 4; i++) {
      kinData[i + 4] = kinData[i];
    }
    // isotropic hardening - turned off by default
    isoData[0] = 0.0;  // b_i
    isoData[1] = 1.0;  // rho_i
    isoData[2] = 0.0;  // b_l
    isoData[3] = 50.0; // R_i
    isoData[4] = 0.0;  // l_yp
    for (unsigned int i = 0; i < 4; i++) {
      isoData[i + 5] = isoData[i];
    }
    // ultimate strength limit - turned off by default (i.e. extremely high
    // limit is applied)
    ultData[0] = 100000000.0 * basicData[0]; // f_u
    ultData[1] = 50.0;                       // R_u
    for (unsigned int i = 0; i < 2; i++) {
      ultData[i + 2] = ultData[i];
    }
    // load history memory - turned on by default
    memData[0] = 50; // cycNum
    // initial stress - zero by default
    initData[0] = 0.0; // sig_0
  }

  // Read the optional string tags to see what else is there to model
  argc = OPS_GetNumRemainingInputArgs();
  // const char *argvLoc = 0;
  while (argc > 1) {
    // char argvLoc[10];
    // if there are no valid string tags throw an error
    const char *argvLoc = OPS_GetString();

    // if the material is asymmetric, modify the number of required parameters
    if (strcmp(argvLoc, "-asym") == 0) {
      numKin = 8;
      numIso = 9;
      numUlt = 4;
    }
    // kinematic hardening
    else if (strcmp(argvLoc, "-kin") == 0) {
      // check if the right number of parameters are provided
      if (OPS_GetDouble(&numKin, kinData) != 0) {
        opserr << "WARNING invalid -kin args for Steel4 (tag: " << iData[0]
               << ")\n"
               << endln;
        return 0;
      }
      // define the unspecified values in case of symmetric hardening
      if (numKin == 4) {
        for (unsigned int i = 0; i < 4; i++) {
          kinData[i + 4] = kinData[i];
        }
      }
    }
    // isotropic hardening with MP characteristics
    else if (strcmp(argvLoc, "-iso") == 0) {
      // check if the right number of parameters are provided
      if (OPS_GetDouble(&numIso, isoData) != 0) {
        opserr << "WARNING invalid -iso args for Steel4 (tag: " << iData[0]
               << ")\n"
               << endln;
        return 0;
      }
      // define the unspecified values in case of symmetric hardening
      if (numIso == 5) {
        for (unsigned int i = 0; i < 4; i++) {
          isoData[i + 5] = isoData[i];
        }
      }
    }
    // ultimate strength limit
    else if (strcmp(argvLoc, "-ult") == 0) {
      // check if the right number of parameters are provided
      if (OPS_GetDouble(&numUlt, ultData) != 0) {
        opserr << "WARNING invalid -ult args for Steel4 (tag: " << iData[0]
               << ")\n"
               << endln;
        return 0;
      }
      // define the unspecified values in case of symmetric hardening
      if (numUlt == 2) {
        for (unsigned int i = 0; i < 2; i++) {
          ultData[i + 2] = ultData[i];
        }
      }
    }
    // load history memory
    else if (strcmp(argvLoc, "-mem") == 0) {
      // check if the right number of parameters are provided
      if (OPS_GetInt(&numMem, memData) != 0) {
        opserr << "WARNING invalid -mem args for Steel4 (tag: " << iData[0]
               << ")\n"
               << endln;
        return 0;
      }
    }
    // initial stress
    else if (strcmp(argvLoc, "-init") == 0) {
      // check if the right number of parameters are provided
      if (OPS_GetDouble(&numInit, initData) != 0) {
        opserr << "WARNING invalid -init args for Steel4 (tag: " << iData[0]
               << ")\n"
               << endln;
        return 0;
      }
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  // Allocate the material
  theMaterial = new Steel4(
      iData[0], basicData[0], basicData[1], kinData[0], kinData[1], kinData[2],
      kinData[3], kinData[4], kinData[5], kinData[6], kinData[7], isoData[0],
      isoData[1], isoData[2], isoData[3], isoData[4], isoData[5], isoData[6],
      isoData[7], isoData[8], ultData[0], ultData[1], ultData[2], ultData[3],
      memData[0], initData[0]);

  // Just in case there was a problem with material creation
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel4\n";
    return 0;
  }

  return theMaterial;
}
