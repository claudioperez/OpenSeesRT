

#include <g3_api.h>


#include <SRC/element/frictionBearing/MultiFP2d.h>
void *OPS_MultiFP2d()
{
  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 3) {
    opserr << "WARNING::MultiFP2d insufficient args\n";
    return theEle;
  }

  // get the id and end nodes
  int iData[5];
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING::MultiFP2d invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];

  numRemainingArgs -= 3;

  const char *nextArg;
  int done = 0;

  int axialCase = 1; //

  opserr << "NUM REMAINING ARGS: " << numRemainingArgs << endln;

  while (done == 0 && numRemainingArgs > 0) {

    nextArg = OPS_GetString();
    //    OPS_GetString(nextArg, 30);
    numRemainingArgs--;

    if (strcmp(nextArg, "-W0") == 0) {
      axialCase = 0;
    } else if (strcmp(nextArg, "-WC") == 0) {
      axialCase = 1;
    } else if (strcmp(nextArg, "-WT") == 0) {
      axialCase = 2;
    }

    else if (strcmp(nextArg, "-material") == 0) {
      if (numRemainingArgs == 3) { // user defined material
        numData = 2;
        if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
          opserr << "WARNING invalid element data\n";
          return 0;
        }

        double dData[1];
        numData = 1;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING error reading element area for element" << eleTag
                 << endln;
          return 0;
        }

        UniaxialMaterial *theFrictionMaterial;
        UniaxialMaterial *theVerticalMaterial;
        theFrictionMaterial = OPS_GetUniaxialMaterial(iData[3]);
        theVerticalMaterial = OPS_GetUniaxialMaterial(iData[4]);

        theEle = new MultiFP2d(eleTag, iData[1], iData[2], theFrictionMaterial,
                               theVerticalMaterial, dData[0], axialCase);
        done = 1;
      } else {
        opserr << "WARNING incorrect #args for MultiFP ele " << eleTag
               << " for -material option" << endln;
      }
    }

    else if (strcmp(nextArg, "-triple") == 0) {
      if (numRemainingArgs == 17) { // triple friction pendulum
        //    R1 R2 R3 h1 h2 h3 D1 D2 D3 d1 d2 d3 mu1 mu2 mu3 Kvert W0
        int type = 3;
        double dData[17];
        numData = 17;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
          opserr << "WARNING error reading element area for element" << eleTag
                 << endln;
          return 0;
        }
        Vector R(3);
        R(0) = dData[0];
        R(1) = dData[1];
        R(2) = dData[2];
        Vector h(3);
        h(0) = dData[3];
        h(1) = dData[4];
        h(2) = dData[5];
        Vector D(3);
        D(0) = dData[6];
        D(1) = dData[7];
        D(2) = dData[8];
        Vector d(3);
        d(0) = dData[9];
        d(1) = dData[10];
        d(2) = dData[11];
        Vector mu(3);
        mu(0) = dData[12];
        mu(1) = dData[13];
        mu(2) = dData[14];

        theEle = new MultiFP2d(eleTag, iData[1], iData[2], type, R, h, D, d, mu,
                               dData[15], dData[16], axialCase);

        done = 1;
      } else {
        opserr << "WARNING incorrect #args for MultiFP ele " << eleTag
               << " for -triple option" << endln;
      }
    } else {
      opserr << "WARNING unknown option: " << nextArg << " for MultiFP ele "
             << eleTag << endln;
      done = 1;
    }

    if (theEle == 0) {
      opserr << "WARNING ran out of memory creating element with tag " << eleTag
             << endln;
      return 0;
    }
  }

  return theEle;
}
