

#include <g3_api.h>


#include <SRC/element/feap/fElmt02.h>
void *OPS_fElmt02()
{

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 2 && ndf != 2) {
    opserr
        << "WARNING - fTruss eleTag? iNode? jNode? A? E? needs ndm=2, ndf=2\n";
    return 0;
  }

  // check the number of arguments is correct
  int argc = OPS_GetNumRemainingInputArgs() + 2;
  int eleArgStart = 2;
  if ((argc - eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element fTruss eleTag? iNode? jNode? A? E?\n";
    return 0;
  }

  // get the id and end nodes
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid truss eleTag, iNode or jNode" << endln;
    return 0;
  }
  int trussId = idata[0], iNode = idata[1], jNode = idata[2];

  double ddata[2];
  numdata = 2;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid truss A or E" << endln;
    return 0;
  }

  double A = ddata[0], E = ddata[1];

  // now create the truss and add it to the Domain
  return new fElmt02(trussId, iNode, jNode, A, E);
}
