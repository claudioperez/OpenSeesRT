

#include <g3_api.h>


#include <SRC/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.h>
void *OPS_ElasticTimoshenkoBeam3d()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs == 0) { // parallel processing
    theElement = new ElasticTimoshenkoBeam3d();
    return theElement;
  }

  if (numRemainingArgs < 11) {
    opserr << "ERROR not enough args provided, want: element "
              "ElasticTimoshenkoBeam3d $tag $iNode $jNode $E $G $A $Jx $Iy $Iz "
              "$Avy $Avz $transTag <-mass $m> <-cMass> \n";
    return 0;
  }

  int numData;
  int iData[5];    // tag, iNode, jNode, transTag, cMass
  double dData[9]; // E, G, A, Jx, Iy, Iz, Avy, Avz, mass

  iData[4] = 0;   // cMass
  dData[8] = 0.0; // mass per unit length

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data (tag, iNode, jNode) element "
              "ElasticTimoshenkoBeam3d.\n";
    return 0;
  }

  numData = 8;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element data (E, G, A, Jx, Iy, Iz, Avy, "
              "Avz) element ElasticTimoshenkoBeam3d "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING invalid element data (transTag) element "
              "ElasticTimoshenkoBeam3d "
           << iData[0] << endln;
    return 0;
  }

  CrdTransf *theTrans = OPS_getCrdTransf(iData[3]);
  if (theTrans == 0) {
    opserr << "WARNING transformation object not found for "
              "ElasticTimoshenkoBeam3d "
           << iData[0] << endln;
    return 0;
  }

  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  while (numRemainingArgs > 0) {
    const char *argvLoc = OPS_GetString();
    numData = 1;

    if ((strcmp(argvLoc, "-mass") == 0) || (strcmp(argvLoc, "mass") == 0) ||
        (strcmp(argvLoc, "-rho") == 0) || (strcmp(argvLoc, "rho") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[8]) != 0) {
        opserr << "WARNING error reading element data (mass) element "
                  "ElasticTimoshenkoBeam3d "
               << iData[0] << endln;
        return 0;
      }
    }
    if ((strcmp(argvLoc, "-lMass") == 0) || (strcmp(argvLoc, "lMass") == 0)) {
      iData[4] = 0; // lumped mass matrix (default)
    }
    if ((strcmp(argvLoc, "-cMass") == 0) || (strcmp(argvLoc, "cMass") == 0)) {
      iData[4] = 1; // consistent mass matrix
    }
    numRemainingArgs = OPS_GetNumRemainingInputArgs();
  }

  theElement = new ElasticTimoshenkoBeam3d(
      iData[0], iData[1], iData[2], dData[0], dData[1], dData[2], dData[3],
      dData[4], dData[5], dData[6], dData[7], *theTrans, dData[8], iData[4]);

  return theElement;
}
