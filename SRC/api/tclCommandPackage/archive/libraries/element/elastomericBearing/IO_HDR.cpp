

#include <g3_api.h>


#include <SRC/element/elastomericBearing/HDR.h>
void *OPS_HDR()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyBearing == 0) {
    opserr << "HDR element - Written by Manish Kumar, University at Buffalo, "
              "2012\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs == 0) { // parallel processing
    theEle = new HDR();
    return theEle;
  }

  if (numArgs != 20 && numArgs != 26 && numArgs != 27 && numArgs != 28 &&
      numArgs != 29 && numArgs != 30 && numArgs != 31 && numArgs != 32) {
    opserr << "ERROR - HDR incorrect # args provided";
    return theEle;
  }

  // get the id and end nodes
  int iData[3];
  double dData[17];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 17;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }

  // get the orientation vector
  Vector x(0);
  Vector y(3);
  y(0) = -1.0;
  y(1) = 0.0;
  y(2) = 0.0;

  // The default values of the parameters
  double kl = 10.0;     // Cavitation parameter
  double phi = 0.5;     // Damage index
  double al = 1.0;      // Strength degradation parameter
  double sDratio = 0.5; // Shear distance ratio
  double m = 0.0;       // Mass of the bearing
  double tc1 = 0.0;     // Cover thickness

  if (numArgs >= 26) {
    double value;
    x.resize(3);
    numData = 1;
    for (int i = 0; i < 3; i++) {
      if (OPS_GetDoubleInput(&numData, &value) != 0) {
        opserr << "WARNING invalid orientation value for element" << eleTag
               << endln;
        return 0;
      } else {
        x(i) = value;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (OPS_GetDoubleInput(&numData, &value) != 0) {
        opserr << "WARNING invalid orientation value for element" << eleTag
               << endln;
        return 0;
      } else {
        y(i) = value;
      }
    }
    if (numArgs >= 27) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &kl) != 0) {
        opserr << "WARNING error reading element property cavitation parameter "
                  "for element"
               << eleTag << endln;
        return 0;
      }
      if (numArgs >= 28) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &phi) != 0) {
          opserr << "WARNING error reading element property damage index for "
                    "element"
                 << eleTag << endln;
          return 0;
        }
        if (numArgs >= 29) {
          numData = 1;
          if (OPS_GetDoubleInput(&numData, &al) != 0) {
            opserr << "WARNING error reading element property strength "
                      "degradation parameter for element"
                   << eleTag << endln;
            return 0;
          }
          if (numArgs >= 30) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
              opserr << "WARNING error reading element property shear distance "
                        "ratio for element"
                     << eleTag << endln;
              return 0;
            }
            if (numArgs >= 31) {
              numData = 1;
              if (OPS_GetDoubleInput(&numData, &m) != 0) {
                opserr
                    << "WARNING error reading element property mass for element"
                    << eleTag << endln;
                return 0;
              }
              if (numArgs == 32) {
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &tc1) != 0) {
                  opserr << "WARNING error reading element property cover "
                            "thickness for element"
                         << eleTag << endln;
                  return 0;
                }
              }
            }
          }
        }
      }
    }
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - HDR \n";
    }
    theEle =
        new HDR(iData[0], iData[1], iData[2], dData[0], dData[1], dData[2],
                dData[3], dData[4], dData[5], int(dData[6]), dData[7], dData[8],
                dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
                dData[15], dData[16], y, x, kl, phi, al, sDratio, m, tc1);
  }

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  return theEle;
}
