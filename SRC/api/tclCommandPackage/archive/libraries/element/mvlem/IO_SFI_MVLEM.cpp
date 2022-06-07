

#include <g3_api.h>


#include <SRC/element/mvlem/SFI_MVLEM.h>
void *OPS_SFI_MVLEM(void)
{
  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  // Parse the script for material parameters
  if (numArgs < 7) {
    opserr << "Want: SFI_MVLEM eleTag Dens iNode jNode m c -thick {fiberThick} "
              "-width {fiberWidth} -rho {Rho} -matConcrete {matTagsConcrete} "
              "-matSteel {matTagsSteel} -matShear {matTagShear}\n";
    return 0;
  }

  int iData[4];
  double dData[1];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid int data for element SFI_MVLEM" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid c for element SFI_MVLEM " << iData[0] << endln;
    return 0;
  }

  int m = iData[3];
  const char *str = 0;

  double *theThickness = new double[m];
  double *theWidth = new double[m];
  int *matTags = new int[m];

  NDMaterial **theMaterials = new NDMaterial *[m];

  numArgs = OPS_GetNumRemainingInputArgs();
  while (numArgs >= (m + 1)) {
    // OPS_GetStringCopy(&str);
    str = OPS_GetString();
    if (strcmp(str, "-thick") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
        opserr << "Invalid thick parameter for SFI_MVLEM   " << iData[0]
               << endln;
        return 0;
      }
    } else if (strcmp(str, "-width") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
        opserr << "Invalid width value for SFI_MVLEM  " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-mat") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
        opserr << "Invalid mat tags for SFI_MVLEM  " << iData[0] << endln;
        return 0;
      }
      for (int i = 0; i < m; i++) {
        theMaterials[i] = 0;
        theMaterials[i] = OPS_getNDMaterial(matTags[i]);
        if (theMaterials[i] == 0) {
          opserr << "Invalid material tag " << matTags[i] << "  for SFI_MVLEM  "
                 << iData[0] << endln;
          return 0;
        }
      }
    }

    // clean up the str
    //    delete [] str;

    numArgs = OPS_GetNumRemainingInputArgs();
  }

  theElement = new SFI_MVLEM(iData[0], iData[1], iData[2], theMaterials,
                             theThickness, theWidth, iData[3], dData[0]);

  // Cleanup dynamic memory
  if (theThickness != 0)
    delete[] theThickness;
  if (theWidth != 0)
    delete[] theWidth;
  if (matTags != 0)
    delete[] matTags;

  if (theMaterials != 0)
    delete[] theMaterials;

  return theElement;
}
