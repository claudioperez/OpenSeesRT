

#include <g3_api.h>


#include <SRC/element/truss/InertiaTruss.h>
OPS_Export void *OPS_InertiaTrussElement()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyTruss == 0) {
    opserr << " \n";
    opserr << "                          InertiaTruss element v1.0\n";
    opserr << "                    by Xiaodong Ji, Yuhao Cheng, Yue Yu\n";
    opserr << "                           Tsinghua University\n";
    opserr << "Please contact jixd@mail.tsinghua.edu.cn, yuhao_cheng@126.com "
              "if anything goes wrong\n";
    opserr << " \n";
    numMyTruss++;
  }

  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr
        << "Invalid Args want: element InertiaTruss $tag $iNode $jNode $mr\n";
    return 0;
  }

  int iData[3];
  double mr = 0.0;
  int ndm = OPS_GetNDM();

  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode) in element "
              "InertiaTruss "
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &mr) != 0) {
    opserr << "WARNING: Invalid mr: element InertiaTruss " << iData[0]
           << " $iNode $jNode $mr\n";
    return 0;
  }

  // now create the InertiaTruss
  theElement = new InertiaTruss(iData[0], ndm, iData[1], iData[2], mr);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element InertiaTruss " << iData[0]
           << " $iNode $jNode $mr\n";
  }

  return theElement;
}
