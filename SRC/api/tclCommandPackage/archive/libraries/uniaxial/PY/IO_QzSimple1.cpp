

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/QzSimple1.h>
void *OPS_QzSimple1(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial QzSimple1 tag? qzType? qult? z50? "
              "suction? c?\n";
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  double ddata[4] = {0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 4)
    numdata = 4;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  UniaxialMaterial *theMaterial = 0;
  theMaterial =
      new QzSimple1(idata[0], idata[1], ddata[0], ddata[1], ddata[2], ddata[3]);

  return theMaterial;
}
