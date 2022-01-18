

#include <g3_api.h>


#include <SRC/material/uniaxial/SLModel.h>
void *OPS_SLModel(G3_Runtime *rt)
{
  // print out some KUDO's
  if (numSLModel == 0) {
    numSLModel++;
    opserr << "SLModel version 2019.2\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[3];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  SLModel tag" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial SLModel tag? Dt?, sgm_ini?, "
              "OP_Material?";
    return 0;
  }

  // create a new material
  theMaterial = new SLModel(iData[0], dData[0], dData[1], dData[2]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SLModel\n";
    return 0;
  }

  // return the material
  return theMaterial;
}
