

#include <g3_api.h>


#include <SRC/element/XMUelements/AV3D4QuadWithSensitivity.h>
void *OPS_AV3D4QuadWithSensitivity(void)
{

  int eleID, numNodes, matTag;
  int nodes[8];
  int i;

  static int idData[6];

  // if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 6) {
    opserr << "element AV3D4Quad incorrect num args .. 6 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element AV3D4Quad error reading integers\n";
    return 0;
  }

  matTag = idData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "command: element AC3D8Hex " << idData[0]
           << " - no NDMaterial with tag " << matTag << " exists\n";
    return 0;
  }

  Element *theEle = new AV3D4QuadWithSensitivity(
      idData[0], idData[1], idData[2], idData[3], idData[4], theMaterial);
  return theEle;
}
