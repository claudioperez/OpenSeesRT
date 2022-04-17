

#include <g3_api.h>


#include <SRC/element/joint/LehighJoint2d.h>
void *OPS_LehighJoint2d()
{
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;

  // check no of arguments
  int numData = OPS_GetNumRemainingInputArgs();
  if (numData != 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LehighJoint eleTag? node1? node2? node3? node4? "
              "matTag1? matTag2? matTag3? ";
    opserr << "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? \n";

    return 0;
  }

  // 1 ele tag, 4 node tags, 9 material tags
  int idata[14];
  int num = 14;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }

  int eleTag = idata[0];
  int ndI = idata[1];
  int ndJ = idata[2];
  int ndK = idata[3];
  int ndL = idata[4];

  UniaxialMaterial *theMats[9];
  for (int i = 0; i < 9; i++) {
    theMats[i] = OPS_getUniaxialMaterial(idata[i + 5]);
    if (theMats[i] == 0) {
      opserr << "WARNING: material not found\n";
      opserr << "Material: " << idata[i + 5];
      opserr << "\nLehighJoint2d element: " << eleTag << endln;
      return 0;
    }
  }

  LehighJoint2d *theEle = 0;
  theEle =
      new LehighJoint2d(eleTag, ndI, ndJ, ndK, ndL, *theMats[0], *theMats[1],
                        *theMats[2], *theMats[3], *theMats[4], *theMats[5],
                        *theMats[6], *theMats[7], *theMats[8]);
  return theEle;
}
