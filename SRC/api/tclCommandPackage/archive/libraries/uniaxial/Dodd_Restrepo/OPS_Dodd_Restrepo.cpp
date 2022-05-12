

#include <g3_api.h>


#include <SRC/material/uniaxial/Dodd_Restrepo.h>
void *OPS_Dodd_Restrepo(void)
{
  if (numDoddRestrepo == 0) {
    numDoddRestrepo++;
    opserr << "Dodd_Restrepo unaxial material - Written by L.L. Dodd & J. "
              "Restepo\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 8 || numArgs > 10) {
    opserr << "WARNING wrong # args: uniaxialMaterial $tag $Fy $Fsu $ESH $ESU "
              "$Youngs $ESHI $FSHI <$OmegaFac>"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[9];
  int numData;

  dData[7] = 1.0; // omegaFac
  dData[8] = 1.0; // Conv

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return 0;
  }

  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;
  }

  theMaterial =
      new Dodd_Restrepo(iData[0], dData[0], dData[1], dData[2], dData[3],
                        dData[4], dData[5], dData[6], dData[7], dData[8]);

  if (theMaterial == 0) {
    opserr
        << "WARNING could not create uniaxialMaterial of type ElasticPPCpp\n";
    return 0;
  }

  return theMaterial;
}
