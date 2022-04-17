

#include <g3_api.h>


#include <SRC/element/HUelements/MultipleNormalSpring.h>
void *OPS_MultipleNormalSpring()
{
  // 3-dim, 6-dof
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleNormalSpring command only works when ndm is 3 "
              "and ndf is 6"
           << endln;
    return 0;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nDivide;

  // arguments (necessary, input with -???)
  int matTag;
  UniaxialMaterial *material;
  int shape;
  double size;

  // arguments (optional, input with -???)
  double lambda = -1.0;
  Vector oriX(0);
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0;
  double mass = 0.0;

  // input comfirmation
  int recvMat = 0;
  int recvShape = 0;
  int recvSize = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass = 0;

  //
  Element *theElement = 0;

  // error flag
  bool ifNoError = true;

  if (OPS_GetNumRemainingInputArgs() <
      4) { // element multipleNormalSpring eleTag? iNode? jNode? nDivide?

    ifNoError = errDetected(ifNoError, "insufficient arguments");

  } else {

    // argv[2~5]
    int idata[4];
    int numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
      ifNoError = errDetected(ifNoError, "invalid int inputs");
    }

    eleTag = idata[0];
    iNode = idata[1];
    jNode = idata[2];
    nDivide = idata[3];
    if (nDivide <= 0) {
      ifNoError = errDetected(ifNoError, "invalid nDivide");
    }

    // argv[6~]
    while (OPS_GetNumRemainingInputArgs() > 0) {

      double value;
      const char *flag = OPS_GetString();
      if (strcmp(flag, "-mat") == 0 &&
          OPS_GetNumRemainingInputArgs() > 0) { // -mat matTag?

        numdata = 1;
        if (OPS_GetIntInput(&numdata, &matTag) < 0) {
          ifNoError = errDetected(ifNoError, "invalid matTag");
        }

        material = OPS_getUniaxialMaterial(matTag);
        if (material == 0) {
          ifNoError = errDetected(ifNoError, "material model not found");
        }

        recvMat++;

      } else if (strcmp(flag, "-shape") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -shape shape?
        const char *shapeflag = OPS_GetString();
        if (strcmp(shapeflag, "round") == 0) {
          shape = 1; // round shape
        } else if (strcmp(shapeflag, "square") == 0) {
          shape = 2; // square
        } else {
          ifNoError = errDetected(
              ifNoError,
              "invalid shape (\"round\" or \"square\" are available)");
        }

        recvShape++;

      } else if (strcmp(flag, "-size") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // -size size?

        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &size) < 0 || size <= 0) {
          ifNoError = errDetected(ifNoError, "invalid size");
        }

        recvSize++;

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
                 OPS_GetNumRemainingInputArgs() > 0) { // <-mass m?> �̓ǂݍ���
        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
          ifNoError = errDetected(ifNoError, "invalid mass");
        }

        recvMass++;

      } else { // invalid option

        ifNoError = errDetected(ifNoError, "invalid optional arguments");
        break;
      }
    }

  } // end input

  // input cofirmation
  // necessary arguments
  if (recvMat != 1) {
    char buf[100];
    sprintf(buf,
            "wrong number of -mat inputs (got %d inputs, but want 1 input)",
            recvMat);
    ifNoError = errDetected(ifNoError, buf);
  }

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

  // optional arguments
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

  // if error detected
  if (!ifNoError) {
    opserr << "------------------------------" << endln;
    // input:
    //  printCommand(argc, argv);
    // want:
    opserr << "Want: element multipleNormalSpring eleTag? iNode? jNode? "
              "nDivide? -mat matTag? -shape shape? -size size? <-lambda "
              "lambda?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    opserr << "========================================" << endln;
    opserr << "" << endln;
    return 0;
  }

  // now create the multipleNormalSpring
  theElement = new MultipleNormalSpring(eleTag, iNode, jNode, nDivide, material,
                                        shape, size, lambda, oriYp, oriX, mass);

  return theElement;
}
