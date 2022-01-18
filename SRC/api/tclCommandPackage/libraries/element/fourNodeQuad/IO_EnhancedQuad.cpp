

#include <g3_api.h>


#include <SRC/element/fourNodeQuad/EnhancedQuad.h>
void *OPS_EnhancedQuad()
{
  if (OPS_GetNDM() != 2 || OPS_GetNDF() != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element ConstantPressureVolumeQuad eleTag? iNode? jNode? "
              "kNode? lNode? thk? type? matTag?\n";
    return 0;
  }

  // EnhancedQuadId, iNode, jNode, kNode, lNode
  int data[5];
  int num = 5;
  if (OPS_GetIntInput(&num, data) < 0) {
    opserr << "WARNING: invalid integer input\n";
    return 0;
  }

  double thk = 1.0;
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
    opserr << "\nConstantPressureVolumeQuad element: " << data[0] << endln;
    return 0;
  }

  return new EnhancedQuad(data[0], data[1], data[2], data[3], data[4], *mat,
                          type, thk);
}
