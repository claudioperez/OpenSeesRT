

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/TzSimple2.h>
void *OPS_TzSimple2(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial TzSimple2 tag? tzType? tult? z50? "
              "dashpot?\n";
    return 0;
  }

  int idata[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  double ddata[3] = {0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 3)
    numdata = 3;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  UniaxialMaterial *theMaterial = 0;
  theMaterial = new TzSimple2(idata[0], MAT_TAG_PySimple1, idata[1], ddata[0],
                              ddata[1], ddata[2]);

  return theMaterial;
}
