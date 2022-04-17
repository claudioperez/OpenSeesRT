

#include <g3_api.h>


#include <SRC/element/elasticBeamColumn/ModElasticBeam2d.h>
void *OPS_ModElasticBeam2d()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numModElasticBeam2d == 0) {
    opserr << "ModElasticBeam2d element -> for Stiffness Modification Factors "
              "by D.Lignos"
           << endln;
    numModElasticBeam2d++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs == 0) { // parallel processing
    theEle = new ModElasticBeam2d();
    return theEle;
  }

  if (numRemainingArgs < 10) {
    opserr << "ERROR not enough args provided, want: element ModElasticBeam2d "
              "tag? iNode? jNode? A? E? I? K11? K33? K44? transfType? <-alpha "
              "$alpha> <-d $d> <-rho $rho> <-cMass>\n";
    return 0;
  }

  int numData;
  int iData[5];    // tag, iNode, jNode, transfTag, cMass
  double dData[9]; // A, E, I, K11, K33, K44, alpha, d, rho

  iData[4] = 0; // cMass
  dData[6] = 0.0;
  dData[7] = 0.0; // alpha and d
  dData[8] = 0.0; // rho

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data (tag, iNode, jNode) element "
              "ElasticBeamColumn2d\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element data (A, E, I) element "
              "ElasticBeamColumn2d "
           << eleTag << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
    opserr << "WARNING error reading element data (K11, K33, K44) element "
              "ElasticBeamColumn2d "
           << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr
        << "WARNING error reading data (transfTag) element ElasticBeamColumn2d "
        << eleTag << endln;
    return 0;
  }

  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  while (numRemainingArgs > 1) {
    const char *argvLoc = OPS_GetString();
    ;
    numData = 1;

    if ((strcmp(argvLoc, "-alpha") == 0) || (strcmp(argvLoc, "-Alpha") == 0) ||
        (strcmp(argvLoc, "-ALPHA") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
        opserr << "WARNING error reading element data (alpha) element "
                  "ElasticBeamColumn2d "
               << eleTag << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-d") == 0) || (strcmp(argvLoc, "-D") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
        opserr << "WARNING error reading element data (D) element "
                  "ElasticBeamColumn2d "
               << eleTag << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-rho") == 0) ||
               (strcmp(argvLoc, "Rho") == 0) ||
               (strcmp(argvLoc, "-RHO") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[8]) != 0) {
        opserr << "WARNING error reading element data (rho) element "
                  "ElasticBeamColumn2d "
               << eleTag << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-lMass") == 0) ||
               (strcmp(argvLoc, "lMass") == 0)) {
      iData[4] = 0; // lumped mass matrix (default)
    } else if ((strcmp(argvLoc, "-cMass") == 0) ||
               (strcmp(argvLoc, "cMass") == 0)) {
      iData[4] = 1; // consistent mass matrix
    }
    numRemainingArgs = OPS_GetNumRemainingInputArgs();
  }

  CrdTransf *theTransf = OPS_getCrdTransf(iData[3]);
  if (theTransf == 0) {
    opserr << "WARNING error could not find a transformation with tag: "
           << iData[3] << "element ElasticBeamColumn2d " << eleTag << endln;
    return 0;
  }

  theEle = new ModElasticBeam2d(
      iData[0], dData[0], dData[1], dData[2], iData[1], iData[2], dData[3],
      dData[4], dData[5], *theTransf, dData[6], dData[7], dData[8], iData[4]);

  //  delete theTransf;

  return theEle;
}
