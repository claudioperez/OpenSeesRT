

#include <g3_api.h>


#include <SRC/element/fourNodeQuad/ConstantPressureVolumeQuad.h>
void *OPS_ConstantPressureVolumeQuad()
{
  if (OPS_GetNDM() != 2 || OPS_GetNDF() != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element ConstantPressureVolumeQuad eleTag? iNode? jNode? "
              "kNode? lNode? thk? matTag?\n";
    return 0;
  }

  // ConstantPressureVolumeQuadId, iNode, jNode, kNode, lNode
  int idata[5];
  int num = 5;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }

  double thk = 1.0;
  num = 1;
  if (OPS_GetDoubleInput(&num, &thk) < 0) {
    opserr << "WARNING: invalid double inputs\n";
    return 0;
  }

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
    opserr << "\nConstantPressureVolumeQuad element: " << idata[0] << endln;
    return 0;
  }

  return new ConstantPressureVolumeQuad(idata[0], idata[1], idata[2], idata[3],
                                        idata[4], *mat, thk);
}
