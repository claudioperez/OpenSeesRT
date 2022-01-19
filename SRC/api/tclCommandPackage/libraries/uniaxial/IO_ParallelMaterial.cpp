

#include <g3_api.h>


#include <SRC/material/uniaxial/ParallelMaterial.h>
void *OPS_ParallelMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  Vector *theFactors = 0;

  int argc = OPS_GetNumRemainingInputArgs();

  if (argc < 2) {
    opserr << "Invalid #args,  want: uniaxialMaterial Parallel $tag $tag1 "
              "$tag2 ... <-factors $fact1 $fact2 ...>"
           << endln;
    return 0;
  }

  // count the number of materials
  int numMats = -1;
  int gotFactors = 0;

  while (argc > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-factors") == 0) {
      gotFactors = 1;
      break;
    }
    numMats++;
    argc = OPS_GetNumRemainingInputArgs();
  }

  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  int numData = 1 + numMats;
  int *iData = new int[numData];
  UniaxialMaterial **theMats = new UniaxialMaterial *[numMats];
  double *dData = 0;

  if (gotFactors) {
    dData = new double[numMats];
    theFactors = new Vector(dData, numMats);
  }

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid data for uniaxialMaterial Parallel" << endln;
    return 0;
  }

  for (int i = 1; i < numMats + 1; i++) {
    UniaxialMaterial *theMat = OPS_getUniaxialMaterial(iData[i]);
    if (theMat == 0) {
      opserr << "WARNING no existing material with tag " << iData[i]
             << " for uniaxialMaterial Parallel" << iData[0] << endln;
      delete[] iData;
      delete[] theMats;
      return 0;
    }
    theMats[i - 1] = theMat;
  }

  if (gotFactors) {
    const char *argvLoc = OPS_GetString();
    if (OPS_GetDoubleInput(&numMats, dData) != 0) {
      opserr << "WARNING invalid factors for uniaxialMaterial Parallel"
             << endln;
      return 0;
    }
  }

  // Parsing was successful, allocate the material
  theMaterial = new ParallelMaterial(iData[0], numMats, theMats, theFactors);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Parallel\n";
    return 0;
  }

  delete[] iData;
  delete[] theMats;

  if (theFactors != 0)
    delete theFactors;

  return theMaterial;
}
