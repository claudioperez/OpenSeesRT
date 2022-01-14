

#include <g3_api.h>


#include <SRC/material/uniaxial/Concrete02.h>
void *OPS_Concrete02(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Concrete02 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 4 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02 " << iData[0]
           << " fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02 " << iData[0]
           << " fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  if (numData == 7)
    theMaterial = new Concrete02(iData[0], dData[0], dData[1], dData[2],
                                 dData[3], dData[4], dData[5], dData[6]);
  else
    theMaterial =
        new Concrete02(iData[0], dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Concrete02 "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
