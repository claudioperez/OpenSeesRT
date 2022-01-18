

#include <g3_api.h>


#include <SRC/material/uniaxial/ElasticBilin.h>
void *OPS_ElasticBilin(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();

  if (argc != 4 && argc != 7) {
    opserr << "WARNING incorrect num args want: uniaxialMaterial ElasticBilin "
              "tag E1P? E2P? eps2P? <E1N? E2N? eps2N?>"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[6];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticBilin tag" << endln;
    return 0;
  }

  argc--;
  numData = argc;
  ;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: uniaxialMaterial ElasticBilin tag "
              "E2P eps2P <E2N? eps2N?>"
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  if (argc == 3)
    theMaterial = new ElasticBilin(iData[0], dData[0], dData[1], dData[2]);
  else
    theMaterial = new ElasticBilin(iData[0], dData[0], dData[1], dData[2],
                                   dData[3], dData[4], dData[5]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type ElasticBilin\n";
    return 0;
  }

  return theMaterial;
}
