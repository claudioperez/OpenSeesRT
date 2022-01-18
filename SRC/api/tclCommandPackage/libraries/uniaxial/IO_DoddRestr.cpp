

#include <g3_api.h>


#include <SRC/material/uniaxial/DoddRestr.h>
void *OPS_DoddRestr(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DoddRestr tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 9 && numData != 12) {
    opserr << "Invalid #args, want: uniaxialMaterial DoddRestr " << iData[0]
           << " Eo? fy? esh? esh1? fsh1? esu? fsu? Pmajor? Pminor? <slcf? "
              "tlcf? Dcrit?>>"
           << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial DoddRestr " << iData[0]
           << " Eo? fy? esh? esh1? fsh1? esu? fsu? Pmajor? Pminor? <slcf? "
              "tlcf? Dcrit?>>"
           << endln;
    return 0;
  }

  if (numData == 9) {
    dData[9] = 0.;
    dData[10] = 0.;
    dData[11] = 0.;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DoddRestr(iData[0], dData[0], dData[1], dData[2], dData[3],
                              dData[4], dData[5], dData[6], dData[7], dData[8],
                              dData[9], dData[10], dData[11]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type DoddRestr "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
