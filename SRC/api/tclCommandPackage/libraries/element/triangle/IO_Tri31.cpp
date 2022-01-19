

#include <g3_api.h>


#include <SRC/element/triangle/Tri31.h>
OPS_Export void *OPS_Tri31()
{
  if (num_Tri31 == 0) {
    num_Tri31++;
    opserr << "Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n";
    // OPS_Error("Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC
    // Berkeley\n",1);
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args, want: element element Tri31 eleTag? iNode? "
              "jNode? kNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return 0;
  }

  int iData[5];
  char *theType;
  double dData[5];
  dData[1] = 0.0;
  dData[2] = 0.0;
  dData[3] = 0.0;
  dData[4] = 0.0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Tri31\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid thickness data: element Tri31 " << iData[0]
           << endln;
    return 0;
  }

  // if (OPS_GetStringCopy(&theType) != 0) {
  //   opserr << "WARNING invalid type, want: ""PlaneStress"" or ""PlaneStrain""
  //   element SSPquad " << iData[0] << endln; return 0;
  // }
  theType = (char *)OPS_GetString();

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer data: element Tri31\n";
    return 0;
  }
  int matID = iData[4];

  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element Tri31 " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  if (numRemainingInputArgs == 11) {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING invalid optional data: element Tri31 " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement =
      new Tri31(iData[0], iData[1], iData[2], iData[3], *theMaterial, theType,
                dData[0], dData[1], dData[2], dData[3], dData[4]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Tri31\n";
    return 0;
  }

  return theElement;
}
