

#include <g3_api.h>


#include <SRC/element/HUelements/MultipleShearSpring.h>
void *OPS_MultipleShearSpring()
{
  // 3-dim, 6-dof
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleShearSpring command only works when ndm is 3 "
              "and ndf is 6"
           << endln;
    return 0;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nSpring;
  int matTag;

  // material
  UniaxialMaterial *material = 0;
  UniaxialMaterial **theMaterials = 0;
  int recvMat = 0;

  // arguments (optional)
  double limDisp = 0.0;
  Vector oriX(0);
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0;
  double mass = 0.0;

  //
  Element *theElement = 0;

  // error flag
  bool ifNoError = true;

  if (OPS_GetNumRemainingInputArgs() <
      6) { // element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat
           // matTag?

    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {

    // argv[2~5]
    int idata[4];
    int numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
      opserr << "WARNING invalid multipleShearSpring int inputs\n";
      ifNoError = false;
    }
    eleTag = idata[0];
    iNode = idata[1];
    jNode = idata[2];
    nSpring = idata[3];
    if (nSpring <= 0) {
      opserr << "WARNING invalid nSpring\n";
      ifNoError = false;
    }

    // argv[6~]
    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char *flag = OPS_GetString();

      double value;

      if (strcmp(flag, "-mat") == 0 &&
          OPS_GetNumRemainingInputArgs() > 0) { // -mat matTag?
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &matTag) < 0) {
          opserr << "WARNING invalid matTag\n";
          ifNoError = false;
        }

        material = OPS_getUniaxialMaterial(matTag);
        if (material == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "multipleShearSpring element: " << eleTag << endln;
          return 0;
        }

        // opserr << "org material " << material->getClassType() << "\n";
        recvMat++;

      } else if (strcmp(flag, "-nMat") == 0 &&
                 OPS_GetNumRemainingInputArgs() >= nSpring) { // -mat matTag?

        theMaterials = new UniaxialMaterial *[nSpring];
        for (int j = 0; j < nSpring; j++) {
          numdata = 1;
          if (OPS_GetIntInput(&numdata, &matTag) < 0) {
            opserr << "WARNING invalid matTag\n";
            ifNoError = false;
          }

          theMaterials[j] = OPS_getUniaxialMaterial(matTag);
          if (theMaterials[j] == 0) {
            opserr << "WARNING material model not found\n";
            opserr << "uniaxialMaterial: " << matTag << endln;
            opserr << "multipleShearSpring element: " << eleTag << endln;
            return 0;
          }
        }
        // opserr << "org material " << material->getClassType() << "\n";
        recvMat++;

      } else if (strcmp(flag, "-orient") == 0 &&
                 OPS_GetNumRemainingInputArgs() >=
                     6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        for (int j = 1; j <= 3; j++) {
          numdata = 1;
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriX(j - 1) = value;
          }
        }

        for (int j = 1; j <= 3; j++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

      } else if (strcmp(flag, "-orient") == 0 &&
                 OPS_GetNumRemainingInputArgs() >=
                     3) { // <-orient yp1? yp2? yp3?> �̓ǂݍ���

        for (int j = 1; j <= 3; j++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

      } else if (strcmp(flag, "-mass") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-mass m?> �̓ǂݍ���
        if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
          opserr << "WARNING invalid mass\n";
          ifNoError = false;
        }

      } else if (strcmp(flag, "-lim") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-lim limDisp?> �̓ǂݍ���
        if (OPS_GetDoubleInput(&numdata, &limDisp) < 0 || limDisp < 0) {
          opserr << "WARNING invalid limDisp\n";
          ifNoError = false;
        }

      } else { // invalid option

        opserr << "WARNING invalid optional arguments \n";
        ifNoError = false;
        break;
      }
    }

  } // end input

  // confirm material
  if (recvMat != 1) {
    opserr << "WARNING wrong number of -mat inputs\n";
    opserr << "got " << recvMat << " inputs, but want 1 input\n";
    ifNoError = false;
  }

  // if error detected
  if (!ifNoError) {
    // input:
    // want:
    opserr << "Want: element multipleShearSpring eleTag? iNode? jNode? "
              "nSpring? -mat matTag? <-lim dsp> <-orient <x1? x2? x3?> yp1? "
              "yp2? yp3?> <-mass m?>\n";
    return 0;
  }

  // now create the multipleShearSpring
  if (theMaterials == 0) {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, nSpring,
                                         material, limDisp, oriYp, oriX, mass);
  } else {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, theMaterials,
                                         nSpring, limDisp, oriYp, oriX, mass);
    delete[] theMaterials;
  }

  return theElement;
}
