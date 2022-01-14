

#include <g3_api.h>


#include <SRC/element/HUelements/YamamotoBiaxialHDR.h>
void *OPS_YamamotoBiaxialHDR()
{
  // 3-dim, 6-dof
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING YamamotoBiaxialHDR command only works when ndm is 3 and "
              "ndf is 6"
           << endln;
    return 0;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  int Tp = 1;
  double DDo;
  double DDi;
  double Hr;

  // arguments (optional)
  double Cr = 1.0;
  double Cs = 1.0;
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

  //
  int numdata = 1;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    // element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?
    // argc =            1           2             3      4      5     6   7 8 9
    // argv =       argv[0]      argv[1]      argv[2]  argv[3] .................
    // argv[8]
    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {

    // argv[2~8]
    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR eleTag\n";
      ifNoError = false;
    }

    // iNode
    if (OPS_GetIntInput(&numdata, &iNode) < 0) {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    // jNode
    if (OPS_GetIntInput(&numdata, &jNode) < 0) {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    // Tp
    const char *tparg = OPS_GetString();
    if (strcmp(tparg, "1") == 0) {
      Tp = 1; // Bridgestone X0.6R (EESD version)
    } else {
      opserr << "WARNING invalid YamamotoBiaxialHDR Tp" << endln;
      ifNoError = false;
    }

    // DDo
    if (OPS_GetDoubleInput(&numdata, &DDo) < 0 || DDo <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDo" << endln;
      ifNoError = false;
    }

    // DDi
    if (OPS_GetDoubleInput(&numdata, &DDi) < 0 || DDi < 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDi" << endln;
      ifNoError = false;
    }

    // Hr
    if (OPS_GetDoubleInput(&numdata, &Hr) < 0 || Hr <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR Hr" << endln;
      ifNoError = false;
    }

    // check print--------------------------------------------/
    //  opserr << "   \n";
    //  opserr << "TclModelBuilder_addYamamotoBiaxialHDR()\n";
    //  opserr << "  tp  = " << Tp << endln;
    //  opserr << "  ddo = " << DDo << endln;
    //  opserr << "  ddi = " << DDi << endln;
    //  opserr << "  hr  = " << Hr << endln;
    //------------------------------------------------------

    // argv[9~]
    while (OPS_GetNumRemainingInputArgs() > 0) {
      double value;
      const char *flag = OPS_GetString();
      if (strcmp(flag, "-orient") == 0 &&
          OPS_GetNumRemainingInputArgs() >=
              6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        // x1, x2, x3
        for (int j = 1; j <= 3; j++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriX(j - 1) = value;
          }
        }

        // yp1, yp2, yp3
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
                     3) { // <-orient yp1? yp2? yp3?>

        for (int j = 1; j <= 3; j++) {
          if (OPS_GetDoubleInput(&numdata, &value) < 0) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

      } else if (strcmp(flag, "-mass") == 0 &&
                 OPS_GetNumRemainingInputArgs() > 0) { // <-mass m?>

        if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
          opserr << "WARNING invalid mass\n";
          ifNoError = false;
        }

      } else if (strcmp(flag, "-coRS") == 0 &&
                 OPS_GetNumRemainingInputArgs() >= 2) { // <-coRS cr? cs?>
        if (OPS_GetDoubleInput(&numdata, &Cr) < 0 || Cr <= 0) {
          opserr << "WARNING invalid cr\n";
          ifNoError = false;
        }
        if (OPS_GetDoubleInput(&numdata, &Cs) < 0 || Cs <= 0) {
          opserr << "WARNING invalid cs\n";
          ifNoError = false;
        }

      } else {

        opserr << "WARNING invalid optional arguments \n";
        ifNoError = false;
        break;
      }
    }

  } // end input

  // if error detected
  if (!ifNoError) {
    // input:
    // want:
    opserr << "Want: element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? "
              "DDi? Hr?  <-coRS cr? cs?> <-orient <x1? x2? x3?> y1? y2? y3?> "
              "<-mass m?>\n";
    return 0;
  }

  // now create the YamamotoBiaxialHDR
  theElement = new YamamotoBiaxialHDR(eleTag, iNode, jNode, Tp, DDo, DDi, Hr,
                                      Cr, Cs, oriYp, oriX, mass);

  // if get here we have successfully created the YamamotoBiaxialHDR and added
  // it to the domain
  return theElement;
}
