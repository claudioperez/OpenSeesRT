

#include <g3_api.h>


#include <SRC/material/uniaxial/ConcretewBeta.h>
void *OPS_ConcretewBeta(void)
{

  int total = OPS_GetNumRemainingInputArgs();

  if (total < 12) {
    opserr << "WARNING incorrect number of arguments\n";
    opserr << "Want: uniaxialMaterial ConcretewBeta $tag $fpc $ec0 $fcint "
              "$ecint $fcres $ecres $ft $ftint $etint $ftres $etres <-lambda "
              "$lambda> <-alpha $alpha> <-beta $bint $ebint $bres $ebres> <-E "
              "$E> <-conf $fcc ecc>\n";
    return 0;
  }

  int tag;
  double dData[11];
  double lambda = 0.5;
  double alpha = 1;
  double bData[4];
  bData[0] = 1;
  bData[1] = 0;
  bData[2] = 1;
  bData[3] = 0;
  double M = 0;
  double E0 = 0;
  double fcc = 0;
  double ecc = 0;

  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
    return 0;
  }

  numData = 11;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
    return 0;
  }

  total -= 12;

  //		char optionFlag[12];
  const char *optionFlag;
  while (total > 0) {
    //			OPS_GetString(optionFlag,12);
    optionFlag = OPS_GetString();
    if (strcmp(optionFlag, "-beta") == 0) {
      numData = 4;
      if (OPS_GetDoubleInput(&numData, bData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of "
                  "-beta for tag "
               << tag << endln;
        return 0;
      }
      total -= 5;
    } else if (strcmp(optionFlag, "-lambda") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &lambda) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of "
                  "-lambda for tag "
               << tag << endln;
        return 0;
      }
      total -= 2;
    } else if (strcmp(optionFlag, "-alpha") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &alpha) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of "
                  "-alpha for tag "
               << tag << endln;
        return 0;
      }
      total -= 2;
    } else if (strcmp(optionFlag, "-M") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &M) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of "
                  "-M for tag "
               << tag << endln;
        return 0;
      }
      total -= 2;
    } else if (strcmp(optionFlag, "-E") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &E0) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of "
                  "-E for tag "
               << tag << endln;
        return 0;
      }
      total -= 2;
    } else if (strcmp(optionFlag, "-conf") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &fcc) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument 1 "
                  "of -conf for tag "
               << tag << endln;
        return 0;
      }
      if (OPS_GetDoubleInput(&numData, &ecc) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument 2 "
                  "of -conf for tag "
               << tag << endln;
        return 0;
      }
      total -= 3;
    } else {
      opserr << "WARNING invalid uniaxialMaterial ConcretewBeta flag " << tag
             << endln;
      return 0;
    }
  }

  UniaxialMaterial *theMaterial = new ConcretewBeta(
      tag, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
      dData[7], dData[8], dData[9], dData[10], lambda, alpha, bData[0],
      bData[1], bData[2], bData[3], M, E0, fcc, ecc);

  return theMaterial;
}
