

#include <g3_api.h>


#include <SRC/material/uniaxial/ConcreteSakaiKawashima.h>
void *OPS_ConcreteSakaiKawashima(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[3];

  int numData = OPS_GetNumRemainingInputArgs();
  if (numData != 4) {
    opserr << "Invalid #args, want: uniaxialMaterial ConcreteSakaiKawashima "
              "E0? sigCC? epsCC?\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteSakaiKawashima tag"
           << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial DoddRestr " << iData[0]
           << " Eo? fy? esh? esh1? fsh1? esu? fsu? Pmajor? Pminor? <slcf? "
              "tlcf? Dcrit?>>"
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new ConcreteSakaiKawashima(iData[0], dData[0], dData[1], dData[2]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ConcreteSakaKawashima  Material\n";
    return 0;
  }

  return theMaterial;
}
