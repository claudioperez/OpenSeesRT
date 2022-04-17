

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthContactNTS2D.h>
void *OPS_ZeroLengthContactNTS2D(void)
{

  if (numZeroLengthContactNTS2D == 0) {
    numZeroLengthContactNTS2D++;
    opserr << "ZeroLengthContactNTS2d - Written by Roozbeh G. Mikola and "
              "N.Sitar, UC Berkeley\n";
  }

  Element *theEle = 0;
  int numData = 0;

  // get the ele tag
  int eleTag, sNdNum, pNdNum;
  numData = 1;

  if (OPS_GetInt(&numData, &eleTag) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied eleTag \n";
    return 0;
  }

  const char *nextString = OPS_GetString();
  if (strcmp(nextString, "-sNdNum") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "
           << "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum "
              "pNode? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  // get the number of secondary nodes
  numData = 1;
  if (OPS_GetInt(&numData, &sNdNum) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied sNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString, "-mNdNum") != 0 &&
      strcmp(nextString, "-pNdNum") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "
           << "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum "
              "pNode? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &pNdNum) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalied sNdNum \n";
    return 0;
  }

  // a quick check on number of args
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3 + sNdNum + pNdNum) {
    opserr << "ZeroLengthContactNTS2D::WARNING too few arguments "
           << "want - element zeroLengthContactNTS2D $tag -sNdNum $sNdNum "
              "-pNdNum $pNdNum -Nodes $Nodes $Kn $Kt $phi";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString, "-Nodes") != 0) {
    opserr << "ZeroLengthContactNTS2D:: expecting "
           << "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum "
              "pNode? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  // read the Nodes values
  numData = sNdNum + pNdNum;
  int *theNodeData = new int[numData];
  ID Nodes(theNodeData, numData);

  if (OPS_GetInt(&numData, theNodeData) != 0) {
    opserr << "ZeroLengthContactNTS2D:: invalid Nodes number value for -Nodes ";
    opserr << eleTag
           << "- element ZeroLengthContactNTS2D eleTag? -sNdNum sNode? -pNdNum "
              "pNode? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  // read the material properties
  numData = 3;
  double dData[3];
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ZeroLengthContactNTS2D::WARNING invalid Kn,Kt or phi\n";
    return 0;
  }

  //
  // now we create the element and add it to the domain
  //

  theEle = new ZeroLengthContactNTS2D(eleTag, sNdNum, pNdNum, Nodes, dData[0],
                                      dData[1], dData[2]);
  return theEle;
}
