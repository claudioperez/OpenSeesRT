

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLength.h>
void *OPS_ZeroLength()
{
  int ndm = OPS_GetNDM();

  //
  // first scan the command line to obtain eleID, iNode, jNode, material ID's
  // and their directions, and the orientation of ele xPrime and yPrime not
  // along the global x and y axis
  //

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 7) {
    opserr << "WARNING too few arguments "
           << "want - element ZeroLength eleTag? iNode? jNode? "
           << "-mat matID1? ... -dir dirMat1? .. "
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";

    return 0;
  }

  // eleTag, iNode, jNode
  int idata[3];
  numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING: failed to get integer data\n";
    return 0;
  }

  // create an array of material pointers, to do this first count
  // the materials to create the array then get matID's and from ModelBuilder
  // obtain pointers to the material objects
  const char *type = OPS_GetString();
  if (strcmp(type, "-mat") != 0) {
    opserr << "WARNING expecting "
           << "- element ZeroLength eleTag? iNode? jNode? "
           << "-mat matID1? ... -dir dirMat1? .. "
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";

    return 0;
  }

  //    std::vector<UniaxialMaterial*> mats;
  // create the array
  ID matTags(0);
  int numMats = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    int mtag;
    numdata = 1;
    // the first one not an int
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (OPS_GetIntInput(&numdata, &mtag) < 0) {
      if (numArgs > OPS_GetNumRemainingInputArgs()) {
        // move current arg back by one
        OPS_ResetCurrentInputArg(-1);
      }
      break;
    }
    matTags[numMats] = mtag;
    numMats++;
  }
  UniaxialMaterial **theMats = new UniaxialMaterial *[numMats];
  UniaxialMaterial **theDampMats = new UniaxialMaterial *[numMats];

  for (int i = 0; i < numMats; i++) {

    theMats[i] = OPS_getUniaxialMaterial(matTags(i));
    theDampMats[i] = 0;

    if (theMats[i] == 0) {
      opserr << "WARNING no material " << matTags(i)
             << "exitsts - element ZeroLength eleTag? iNode? jNode? "
             << "-mat matID1? ... -dir dirMat1? .. "
             << "<-orient x1? x2? x3? y1? y2? y3?>\n";
      return 0;
    }
  }

  // now read the dirn ID's for the materials added
  type = OPS_GetString();
  if (strcmp(type, "-dir") != 0 && strcmp(type, "-dof") != 0) {
    opserr << "WARNING expecting -dir flag "
           << "- element ZeroLength eleTag? iNode? jNode? "
           << "-mat matID1? ... -dir dirMat1? .. "
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < numMats) {
    opserr << "WARNING not enough directions provided for ele " << idata[0]
           << "- element ZeroLength eleTag? iNode? jNode? "
           << "-mat matID1? ... -dir dirMat1? .. "
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";
    return 0;
  }

  ID dirs(numMats);
  if (OPS_GetIntInput(&numMats, &dirs(0)) < 0) {
    opserr << "WARNING invalid dir\n";
    return 0;
  }
  for (int i = 0; i < dirs.Size(); i++) {
    dirs(i)--; // subscrit to C++
  }

  // create the vectors for the element orientation
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;

  // finally check the command line to see if user specified orientation
  int doRayleighDamping = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    type = OPS_GetString();
    if (strcmp(type, "-doRayleigh") == 0) {
      doRayleighDamping = 1;
      if (OPS_GetNumRemainingInputArgs() > 0) {
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &doRayleighDamping) < 0) {
          opserr << "WARNING: invalid integer\n";
          return 0;
        }
      }
    } else if (strcmp(type, "-dampMats") == 0) {
      doRayleighDamping = 2;
      numdata = 1;
      int matType;
      for (int i = 0; i < numMats; i++) {
        // the first one not an int
        if (OPS_GetIntInput(&numdata, &matType) < 0) {
          UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matType);
          if (theMat == 0) {
            opserr << "WARNING no damp material material " << matType
                   << " for zeroLength ele: " << idata[0] << endln;
            return 0;
          } else {
            theDampMats[i] = theMat;
          }
        }
      }

    } else if (strcmp(type, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING: insufficient orient values\n";
        return 0;
      }
      numdata = 3;
      if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
        opserr << "WARNING: invalid double input\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
        opserr << "WARNING: invalid double input\n";
        return 0;
      }
    }
  }

  Element *theEle = 0;
  if (doRayleighDamping != 2)
    theEle = new ZeroLength(idata[0], ndm, idata[1], idata[2], x, y, numMats,
                            theMats, dirs, doRayleighDamping);
  else
    theEle = new ZeroLength(idata[0], ndm, idata[1], idata[2], x, y, numMats,
                            theMats, theDampMats, dirs, doRayleighDamping);

  // return the memory we stole and return OK
  delete[] theMats;
  delete[] theDampMats;

  return theEle;
}
