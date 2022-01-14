

#include <g3_api.h>


#include <SRC/element/mvlem/MVLEM_3D.h>
void *OPS_MVLEM_3D(void)
{
  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  // Parse the script for material parameters
  if (numArgs < 18) {
    opserr
        << "Want: MVLEM_3D eleTag iNode jNode kNode lNode m -thick "
           "{fiberThick} -width {fiberWidth} -rho {Rho} -matConcrete "
           "{matTagsConcrete} -matSteel {matTagsSteel} -matShear {matTagShear} "
           "<-CoR c> <-ThickMod tMod> <-Poisson nu> <-Density Dens>\n";
    return 0;
  }

  int iData[6];
  double dData[4];

  // set defaults
  dData[0] = 0.4; // c
  dData[1] =
      0.63; // tMod (equivalent to cracked out-of-plan stiffness of 0.25Ig)
  dData[2] = 0.25; // Poisson (concrete)
  dData[3] = 0.0;  // Density

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for element MVLEM_3D" << endln;
    return 0;
  }

  numData = 5;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING iNode jNode kNode lNode or m for element MVLEM_3D"
           << iData[0] << endln;
    return 0;
  }

  int m = iData[5];
  const char *str = 0;

  double *theThickness = new double[m];
  double *theWidth = new double[m];
  double *theRho = new double[m];
  int *matTags = new int[m];

  UniaxialMaterial **theMaterialsConcrete = new UniaxialMaterial *[m];
  UniaxialMaterial **theMaterialsSteel = new UniaxialMaterial *[m];
  UniaxialMaterial **theMaterialsShear = new UniaxialMaterial *[1];

  numArgs = OPS_GetNumRemainingInputArgs();
  while (numArgs > 0) {
    // OPS_GetStringCopy(&str);
    str = OPS_GetString();
    if (strcmp(str, "-thick") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
        opserr << "Invalid thick parameter for MVLEM   " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-width") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
        opserr << "Invalid width value for MVLEM  " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-rho") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theRho) != 0) {
        opserr << "Invalid width value for MVLEM  " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-matConcrete") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
        opserr << "Invalid width value for MVLEM  " << iData[0] << endln;
        return 0;
      }
      for (int i = 0; i < m; i++) {
        theMaterialsConcrete[i] = 0;
        theMaterialsConcrete[i] = OPS_getUniaxialMaterial(matTags[i]);
        if (theMaterialsConcrete[i] == 0) {
          opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  "
                 << iData[0] << endln;
          return 0;
        }
      }
    } else if (strcmp(str, "-matSteel") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
        opserr << "Invalid steel tags for MVLEM  " << iData[0] << endln;
        return 0;
      }
      for (int i = 0; i < m; i++) {
        theMaterialsSteel[i] = 0;
        theMaterialsSteel[i] = OPS_getUniaxialMaterial(matTags[i]);
        if (theMaterialsSteel[i] == 0) {
          opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  "
                 << iData[0] << endln;
          return 0;
        }
      }
    } else if (strcmp(str, "-matShear") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
        opserr << "Invalid shear tags for MVLEM  " << iData[0] << endln;
        return 0;
      }
      for (int i = 0; i < 1; i++) {
        theMaterialsShear[i] = 0;
        theMaterialsShear[i] = OPS_getUniaxialMaterial(matTags[i]);
        if (theMaterialsShear[i] == 0) {
          opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  "
                 << iData[0] << endln;
          return 0;
        }
      }
    }
    // optional parameters
    else if (strcmp(str, "-CoR") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
        opserr << "Invalid CoR parameter for MVLEM   " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-ThickMod") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
        opserr << "Invalid thickMod parameter for MVLEM   " << iData[0]
               << endln;
        return 0;
      }
    } else if (strcmp(str, "-Poisson") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
        opserr << "Invalid Poisson parameter for MVLEM   " << iData[0] << endln;
        return 0;
      }
    } else if (strcmp(str, "-Density") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
        opserr << "Invalid Dens parameter for MVLEM   " << iData[0] << endln;
        return 0;
      }
    }
    numArgs = OPS_GetNumRemainingInputArgs();
  }

  theElement = new MVLEM_3D(iData[0], dData[3], iData[1], iData[2], iData[3],
                            iData[4], theMaterialsConcrete, theMaterialsSteel,
                            theMaterialsShear, theRho, theThickness, theWidth,
                            iData[5], dData[0], dData[2], dData[1]);

  // Cleanup dynamic memory
  if (theThickness != 0)
    delete[] theThickness;
  if (theWidth != 0)
    delete[] theWidth;
  if (theRho != 0)
    delete[] theRho;
  if (matTags != 0)
    delete[] matTags;

  if (theMaterialsConcrete != 0)
    delete[] theMaterialsConcrete;
  if (theMaterialsSteel != 0)
    delete[] theMaterialsSteel;
  if (theMaterialsShear != 0)
    delete[] theMaterialsShear;

  return theElement;
}
