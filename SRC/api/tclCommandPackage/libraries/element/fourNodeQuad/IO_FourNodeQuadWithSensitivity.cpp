

#include <g3_api.h>


#include <SRC/element/fourNodeQuad/FourNodeQuadWithSensitivity.h>
void *OPS_FourNodeQuadWithSensitivity()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm != 2 || ndf != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element FourNodeQuadWithSensitivity eleTag? iNode? jNode? "
              "kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return 0;
  }

  // FourNodeQuadId, iNode, jNode, kNode, lNode
  int idata[5];
  int num = 5;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  double thk = 0.0;
  num = 1;
  if (OPS_GetDoubleInput(&num, &thk) < 0) {
    opserr << "WARNING: invalid double inputs\n";
    return 0;
  }

  const char *type = OPS_GetString();

  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING: invalid matTag\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << matTag;
    opserr << "\nFourNodeQuad element: " << idata[0] << endln;
    return 0;
  }

  // p, rho, b1, b2
  double data[4] = {0, 0, 0, 0};
  num = OPS_GetNumRemainingInputArgs();
  if (num > 4) {
    num = 4;
  }
  if (num > 0) {
    if (OPS_GetDoubleInput(&num, data) < 0) {
      opserr << "WARNING: invalid integer data\n";
      return 0;
    }
  }

  return new FourNodeQuadWithSensitivity(idata[0], idata[1], idata[2], idata[3],
                                         idata[4], *mat, type, thk, data[0],
                                         data[1], data[2], data[3]);
}
