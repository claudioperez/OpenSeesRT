

#include <g3_api.h>


#include <SRC/element/fourNodeQuad/FourNodeQuad3d.h>
void *OPS_FourNodeQuad3d()
{

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new FourNodeQuad3d();
    return theEle;
  }

  if (numRemainingArgs != 8 && numRemainingArgs != 12) {
    opserr << "ERROR - FourNodeQuad3d not enough args provided, want: element "
              "FourNodeQuad3d tag? iNode? jNode? kNode? lNode? thickness? "
              "type? matID? <p? rho? b1? b2?>\n";
  }

  // get the id and end nodes
  int iData[6];
  double dData[5];
  dData[1] = 0.0;
  dData[2] = 0.0;
  dData[3] = 0.0;
  dData[4] = 0.0;

  int numData;
  int matTag = 0;
  int eleTag = 0;
  const char *pType;

  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid element data\n";
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid thickness for element: "
           << eleTag << "\n";
    return 0;
  }

  pType = OPS_GetString();
  if (pType != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid pType for element: "
           << eleTag << "\n";
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid matTag for element: "
           << eleTag << "\n";
    delete[] pType;
    return 0;
  }

  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element "
           << eleTag << endln;
    return 0;
  }

  if (numRemainingArgs == 12) {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING element FourNodeQuad3d : invalid optional args for "
                "element: "
             << eleTag << "\n";
      delete[] pType;
      return 0;
    }
  }

  // now create the truss and add it to the Domain
  theEle = new FourNodeQuad3d(eleTag, iData[1], iData[2], iData[3], iData[4],
                              *theMaterial, pType, dData[0], dData[1], dData[2],
                              dData[3], dData[4]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    delete theMaterial;
    delete[] pType;
    return 0;
  }

  delete[] pType;
  return theEle;
}
