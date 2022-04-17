

#include <g3_api.h>


#include <SRC/material/uniaxial/PY/PySimple3.h>
void *OPS_PySimple3(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial PySimple3 tag? pult? pyield? ke? C? "
              "dashpot? "
           << endln;
    return 0;
  }

  int idata[1];
  numdata = 1;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }

  double ddata[5] = {0, 0, 0, 0, 0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 5)
    numdata = 5;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }

  UniaxialMaterial *theMaterial = 0;
  theMaterial = new PySimple3(idata[0], MAT_TAG_PySimple3, ddata[0], ddata[1],
                              ddata[2], ddata[3], ddata[4]);

  return theMaterial;
}
