

#include <g3_api.h>


#include <SRC/element/elastomericBearing/ElastomericX.h>
void *OPS_ElastomericX()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyBearing == 0) {
    opserr << "ElastomericX element - Written by Manish Kumar, University at "
              "Buffalo, 2012\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs == 0) { // parallel processing
    theEle = new ElastomericX();
    return theEle;
  }

  if (numArgs != 12 && numArgs != 18 && numArgs != 19 && numArgs != 20 &&
      numArgs != 21 && numArgs != 22 && numArgs != 23 && numArgs != 24 &&
      numArgs != 25 && numArgs != 26 && numArgs != 27 && numArgs != 28 &&
      numArgs != 29) {
    opserr << "ERROR - ElastomericX incorrect # args provided";
    return theEle;
  }

  // get the id and end nodes
  int iData[3];
  double dData[9];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 9;
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

  // get the tags of the properties
  int tag1 = 0; // Cavitation and post-cavitation
  int tag2 = 0; // Buckling load variation
  int tag3 = 0; // Horizontal stiffness variation
  int tag4 = 0; // Vertical stiffness variation

  // The default values of the parameters
  double kl = 10.0;     // Cavitation parameter
  double phi = 0.5;     // Damage index
  double al = 1.0;      // Strength degradation parameter
  double sDratio = 0.5; // Shear distance ratio
  double m = 0.0;       // Mass of the bearing
  double cd1 = 0.0;     // Viscous damping parameter
  double tc1 = 0.0;     // Cover thickness

  if (numArgs >= 18) {
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
    if (numArgs >= 19) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &kl) != 0) {
        opserr << "WARNING error reading element property cavitation parameter "
                  "for element"
               << eleTag << endln;
        return 0;
      }
      if (numArgs >= 20) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &phi) != 0) {
          opserr << "WARNING error reading element property damage index for "
                    "element"
                 << eleTag << endln;
          return 0;
        }
        if (numArgs >= 21) {
          numData = 1;
          if (OPS_GetDoubleInput(&numData, &al) != 0) {
            opserr << "WARNING error reading element property strength "
                      "degradation parameter for element"
                   << eleTag << endln;
            return 0;
          }
          if (numArgs >= 22) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
              opserr << "WARNING error reading element property shear distance "
                        "ratio for element"
                     << eleTag << endln;
              return 0;
            }
            if (numArgs >= 23) {
              numData = 1;
              if (OPS_GetDoubleInput(&numData, &m) != 0) {
                opserr
                    << "WARNING error reading element property mass for element"
                    << eleTag << endln;
                return 0;
              }
              if (numArgs >= 24) {
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &cd1) != 0) {
                  opserr << "WARNING error reading element property viscous "
                            "damping parameter for element"
                         << eleTag << endln;
                  return 0;
                }
                if (numArgs >= 25) {
                  numData = 1;
                  if (OPS_GetDoubleInput(&numData, &tc1) != 0) {
                    opserr << "WARNING error reading element property cover "
                              "thickness for element"
                           << eleTag << endln;
                    return 0;
                  }
                  if (numArgs >= 26) {
                    numData = 1;
                    if (OPS_GetIntInput(&numData, &tag1) != 0) {
                      opserr << "WARNING error reading element properties for "
                                "element"
                             << eleTag << endln;
                      return 0;
                    }
                    if (numArgs >= 27) {
                      numData = 1;
                      if (OPS_GetIntInput(&numData, &tag2) != 0) {
                        opserr << "WARNING error reading element properties "
                                  "for element"
                               << eleTag << endln;
                        return 0;
                      }
                      if (numArgs >= 28) {
                        numData = 1;
                        if (OPS_GetIntInput(&numData, &tag3) != 0) {
                          opserr << "WARNING error reading element properties "
                                    "for element"
                                 << eleTag << endln;
                          return 0;
                        }
                        if (numArgs == 29) {
                          numData = 1;
                          if (OPS_GetIntInput(&numData, &tag4) != 0) {
                            opserr << "WARNING error reading element "
                                      "properties for element"
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
      opserr << ", for space problem need 6 - ElastomericX \n";
    }
    theEle = new ElastomericX(iData[0], iData[1], iData[2], dData[0], dData[1],
                              dData[2], dData[3], dData[4], dData[5], dData[6],
                              dData[7], dData[8], y, x, kl, phi, al, sDratio, m,
                              cd1, tc1, tag1, tag2, tag3, tag4);
  }

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    return 0;
  }

  return theEle;
}
