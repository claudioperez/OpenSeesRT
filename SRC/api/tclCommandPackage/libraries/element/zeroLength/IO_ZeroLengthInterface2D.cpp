

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthInterface2D.h>
void *OPS_ZeroLengthInterface2D(void)
{

  if (numZeroLengthInterface2D == 0) {
    numZeroLengthInterface2D++;
    opserr << "ZeroLengthContactNTS2d - Written by Roozbeh G. Mikola and "
              "N.Sitar, UC Berkeley\n";
  }

  Element *theEle = 0;
  int numData = 0;

  // get the ele tag
  int eleTag, sNdNum, pNdNum, sDOF, mDOF;
  numData = 1;

  if (OPS_GetInt(&numData, &eleTag) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalid eleTag \n";
    return 0;
  }

  const char *nextString = OPS_GetString();

  if (strcmp(nextString, "-sNdNum") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -sNdNum \n";
    return 0;
  }

  // get the number of secondary nodes
  numData = 1;
  if (OPS_GetInt(&numData, &sNdNum) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied sNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString, "-mNdNum") != 0 &&
      strcmp(nextString, "-pNdNum") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -pNdNum\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &pNdNum) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied pNdNum \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString, "-dof") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -sdof in "
           << "element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum "
              "pNdNum? -dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &sDOF) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied sDOF\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &mDOF) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalied mDOF\n";
    return 0;
  }

  // a quick check on number of args
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3 + sNdNum + pNdNum) {
    opserr << "ZeroLengthInterface2D::WARNING too few arguments "
           << "element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum "
              "pNdNum? -dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi? \n";
    return 0;
  }

  numData = 10;
  nextString = OPS_GetString();

  if (strcmp(nextString, "-Nodes") != 0) {
    opserr << "ZeroLengthInterface2D:: expecting -Nodes\n";
    return 0;
  }

  // read the Nodes values
  numData = sNdNum + pNdNum;
  int *theNodeData = new int[numData];
  ID Nodes(theNodeData, numData);

  if (OPS_GetInt(&numData, theNodeData) != 0) {
    opserr << "ZeroLengthInterface2D:: not enough node tags provided for ele: ";
    opserr << eleTag << "\n";
    return 0;
  }

  // read the material properties
  numData = 3;
  double dData[3];
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ZeroLengthInterface2D::WARNING invalid Kn,Kt or phi\n";
    return 0;
  }

  //
  // now we create the element and add it to the domain
  //

  theEle = new ZeroLengthInterface2D(eleTag, sNdNum, pNdNum, sDOF, mDOF, Nodes,
                                     dData[0], dData[1], dData[2]);
  return theEle;
}
