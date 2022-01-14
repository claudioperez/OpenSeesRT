

#include <g3_api.h>


#include <SRC/element/XMUelements/AC3D8HexWithSensitivity.h>
void *OPS_AC3D8HexWithSensitivity(void)
{

  int matTag;

  static int idData[10];

  // if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 10) {
    opserr << "element AC3D8Hex incorrect num args .. 11 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element AC3D8Hex error reading integers\n";
    return 0;
  }

  matTag = idData[9];
  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "command: element AC3D8Hex " << idData[0]
           << " - no NDMaterial with tag " << matTag << " exists\n";
    return 0;
  }

  Element *theEle = new AC3D8HexWithSensitivity(
      idData[0], idData[1], idData[2], idData[3], idData[4], idData[5],
      idData[6], idData[7], idData[8], theMaterial);
  return theEle;
}
