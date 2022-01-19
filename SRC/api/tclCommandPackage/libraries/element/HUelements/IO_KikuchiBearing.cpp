

#include <g3_api.h>


#include <SRC/element/HUelements/KikuchiBearing.h>
void *OPS_KikuchiBearing()
{

  // 3-dim, 6dof
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING KikuchiBearing command only works when ndm is 3 and ndf "
              "is 6"
           << endln;
    return 0;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  // arguments (necessary, input with -???)
  int shape;
  double size;
  double totalRubber;
  int nMSS;
  int matMSSTag;
  UniaxialMaterial *matMSS;
  int nMNS;
  int matMNSTag;
  UniaxialMaterial *matMNS;

  // arguments (optional, input with -???)
  double totalHeight = -1.0; // default: Norm(I->J)
  double limDisp = -1.0;     // default: INF
  double lambda = -1.0;      // default: INF
  Vector oriX(0);            // default: local-x Vec(I->J)
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0; // default: global-Y
  double mass = 0.0;
  bool ifPDInput = true;
  bool ifTilt = true;
  double adjCi = 0.5;
  double adjCj = 0.5;
  bool ifBalance = false;
  double limFo = -1.0; // default: INF
  double limFi = -1.0; // default: INF
  int nIter = 1;

  // input comfirmation
  int recvShape = 0;
  int recvSize = 0;
  int recvHeight = 0;
  int recvNMSS = 0;
  int recvMatMSS = 0;
  int recvLimDisp = 0;
  int recvNMNS = 0;
  int recvMatMNS = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass = 0;
  int recvIfPD = 0;
  int recvIfTl = 0;
  int recvAdj = 0;
  int recvBal = 0;

  //
  Element *theElement = 0;

  // error flag
  bool ifNoError = true;

  if (OPS_GetNumRemainingInputArgs() <
      3) { // element KikuchiBearing eleTag? iNode? jNode?

    ifNoError = errDetected(ifNoError, "insufficient arguments");

  } else {

    // argv[2~4]
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
      ifNoError = errDetected(ifNoError, "invalid eleTag");
    }

    if (OPS_GetIntInput(&numdata, &iNode) < 0) {
      ifNoError = errDetected(ifNoError, "invalid iNode");
    }

    if (OPS_GetIntInput(&numdata, &jNode) < 0) {
      ifNoError = errDetected(ifNoError, "invalid jNode");
    }

    // argv[5~]
    while (OPS_GetNumRemainingInputArgs() > 0) {

      double value;
      const char *flag = OPS_GetString();

      if (strcmp(flag, "-shape") == 0 &&
          OPS_GetNumRemainingInputArgs() > 0) { // -shape shape?
        const char *shapeflag = OPS_GetString();
        if (strcmp(shapeflag, "round") == 0) {
          shape = 1; // round
        } else if (strcmp(shapeflag, "square") == 0) {
          shape = 2; // square
        } else {
          ifNoError = errDetected(
              ifNoError,
              "invalid shape (\"round\" or \"square\" are available)");
        }

        recvShape++;

      } else if (strcmp(flag, "-size") == 0 &&
                 OPS_GetNumRemainingInputArgs() >
                     1) { // -size size? totalRubber?

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &size) < 0 || size <= 0) {
          ifNoError = errDetected(ifNoError, "invalid size");
        }

        if (OPS_GetDoubleInput(&numdata, &totalRubber) < 0 ||
            totalRubber <= 0) {
          ifNoError = errDetected(ifNoError, "invalid totalRubber");
        }

        recvSize++;

      } else if (strcmp(flag, "-totalHeight") == 0 &&
                 OPS_GetNumRemainingInputArgs() >
                     0) { // -totalHeight totalHeight?

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &totalHeight) < 0 ||
            totalHeight <= 0) {
          ifNoError = errDetected(ifNoError, "invalid totalHeight");
        }

        recvHeight++;

      } else if (strcmp(flag, "-nMSS") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -nMSS nMSS?

        numdata = 1;
        if (OPS_GetIntInput(&numdata, &nMSS) < 0 || nMSS <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nMSS");
        }

        recvNMSS++;

      } else if (strcmp(flag, "-matMSS") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -matMSS matMSSTag?

        numdata = 1;
        if (OPS_GetIntInput(&numdata, &matMSSTag) < 0) {
          ifNoError = errDetected(ifNoError, "invalid matMSSTag");
        }

        matMSS = OPS_getUniaxialMaterial(matMSSTag);
        if (matMSS == 0) {
          ifNoError =
              errDetected(ifNoError, "material for MSS model not found");
        }

        recvMatMSS++;

      } else if (strcmp(flag, "-limDisp") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-limDisp limDisp?>

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &limDisp) < 0 || limDisp < 0) {
          ifNoError = errDetected(ifNoError, "invalid limDisp");
        }

        recvLimDisp++;

      } else if (strcmp(flag, "-nMNS") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -nMNS nMNS?

        numdata = 1;
        if (OPS_GetIntInput(&numdata, &nMNS) < 0 || nMNS <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nMNS");
        }

        recvNMNS++;

      } else if (strcmp(flag, "-matMNS") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -matMNS matMNSTag?
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &matMNSTag) < 0) {
          ifNoError = errDetected(ifNoError, "invalid matMNSTag");
        }

        matMNS = OPS_getUniaxialMaterial(matMNSTag);
        if (matMNS == 0) {
          ifNoError =
              errDetected(ifNoError, "material for MNS model not found");
        }

        recvMatMNS++;

      } else if (strcmp(flag, "-lambda") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-lambda lambda?>
        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &lambda) < 0 || lambda < 0) {
          ifNoError = errDetected(ifNoError, "invalid lambda");
        }

        recvLambda++;

      } else if (strcmp(flag, "-orient") == 0 &&
                 OPS_GetNumRemainingInputArgs() >=
                     6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        for (int j = 1; j <= 3; j++) {
          numdata = 1;
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriX(j - 1) = value;
          }
        }

        for (int j = 1; j <= 3; j++) {
          numdata = 1;
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;

      } else if (strcmp(flag, "-orient") == 0 &&
                 OPS_GetNumRemainingInputArgs() >=
                     3) { // <-orient yp1? yp2? yp3?>

        for (int j = 1; j <= 3; j++) {
          numdata = 1;
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;

      } else if (strcmp(flag, "-mass") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-mass mass?>

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
          ifNoError = errDetected(ifNoError, "invalid mass");
        }

        recvMass++;

      } else if (strcmp(flag, "-noPDInput") == 0) { // <-noPDInput>

        ifPDInput = false;

        recvIfPD++;

      } else if (strcmp(flag, "-noTilt") == 0) { // <-noTilt>

        ifTilt = false;

        recvIfTl++;

      } else if (strcmp(flag, "-adjustPDOutput") == 0 &&
                 OPS_GetNumRemainingInputArgs() >
                     1) { // -adjustPDOutput ci? cj?

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &adjCi) < 0 || adjCi <= 0) {
          ifNoError = errDetected(ifNoError, "invalid ci");
        }

        if (OPS_GetDoubleInput(&numdata, &adjCj) < 0 || adjCj <= 0) {
          ifNoError = errDetected(ifNoError, "invalid cj");
        }

        recvAdj++;

      } else if (strcmp(flag, "-doBalance") == 0 &&
                 OPS_GetNumRemainingInputArgs() >
                     2) { // -doBalance limFo? limFi? nIter?

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &limFo) < 0 || limFo <= 0) {
          ifNoError = errDetected(ifNoError, "invalid limFo");
        }

        if (OPS_GetDoubleInput(&numdata, &limFi) < 0 || limFi <= 0) {
          ifNoError = errDetected(ifNoError, "invalid limFi");
        }

        if (OPS_GetIntInput(&numdata, &nIter) < 0 || nIter <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nIter");
        }

        ifBalance = true;

        recvBal++;

      } else { // invalid option

        ifNoError = errDetected(ifNoError, "invalid optional arguments");
        break;
      }
    }

  } // end input

  // input cofirmation
  // necessary arguments
  if (recvShape != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -shape inputs (got %d inputs, but want 1 input)",
            recvShape);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvSize != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -size inputs (got %d inputs, but want 1 input)",
            recvSize);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvNMSS != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -NMSS inputs (got %d inputs, but want 1 input)",
            recvNMSS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMatMSS != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -matMSS inputs (got %d inputs, but want 1 input)",
            recvMatMSS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvNMNS != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -NMNS inputs (got %d inputs, but want 1 input)",
            recvNMNS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMatMNS != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -matMNS inputs (got %d inputs, but want 1 input)",
            recvMatMNS);
    ifNoError = errDetected(ifNoError, buf);
  }

  // optional arguments
  if (recvHeight >= 2) {
    char buf[100];
    sprintf(
        buf,
        "wrong number of -totalHeight inputs (got %d inputs, but want 1 input)",
        recvHeight);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvLimDisp >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -limDisp inputs (got %d inputs, but want 1 input)",
            recvLimDisp);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvLambda >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -lambda inputs (got %d inputs, but want 1 input)",
            recvLambda);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvOrient >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -ori inputs (got %d inputs, but want 1 input)",
            recvOrient);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMass >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -mass inputs (got %d inputs, but want 1 input)",
            recvMass);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvIfPD >= 2) {
    char buf[100];
    sprintf(
        buf,
        "wrong number of -noPDInput inputs (got %d inputs, but want 1 input)",
        recvIfPD);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvIfTl >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -noTilt inputs (got %d inputs, but want 1 input)",
            recvIfTl);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvAdj >= 2) {
    char buf[100];
    sprintf(buf,
            "wrong number of -adjustPDOutput inputs (got %d inputs, but want 1 "
            "input)",
            recvAdj);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvBal >= 2) {
    char buf[100];
    sprintf(
        buf,
        "wrong number of -doBalance inputs (got %d inputs, but want 1 input)",
        recvBal);
    ifNoError = errDetected(ifNoError, buf);
  }

  // if error detected
  if (!ifNoError) {
    opserr << "------------------------------" << endln;
    // input:
    // printCommand(argc, argv);
    // want:
    opserr << "Want: element KikuchiBearing eleTag? iNode? jNode?\n";
    opserr << "                             -shape shape? -size size? "
              "totalRubber? <-totalHeight totalHeight?>\n";
    opserr << "                             -nMSS nMSS? -matMSS matMSSTag? "
              "<-lim limDisp?>\n";
    opserr << "                             -nMNS nMNS? -matMNS matMNSTag? "
              "<-lambda lambda?>\n";
    opserr << "                             <-orient <x1? x2? x3?> yp1? yp2? "
              "yp3?> <-mass m?>\n";
    opserr << "                             <-noPDInput> <-noTilt> "
              "<-adjustPDOutput ci? cj?> <-doBalance limFo? limFi? nIter?>\n";
    opserr << "========================================" << endln;
    opserr << "" << endln;
    return 0;
  }

  // now create the KikuchiBearing
  theElement = new KikuchiBearing(
      eleTag, iNode, jNode, shape, size, totalRubber, totalHeight, nMSS, matMSS,
      limDisp, nMNS, matMNS, lambda, oriYp, oriX, mass, ifPDInput, ifTilt,
      adjCi, adjCj, ifBalance, limFo, limFi, nIter);
  return theElement;
}
