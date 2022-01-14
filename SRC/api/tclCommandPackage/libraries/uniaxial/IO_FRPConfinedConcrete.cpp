

#include <g3_api.h>


#include <SRC/material/uniaxial/FRPConfinedConcrete.h>
void *OPS_FRPConfinedConcrete(void)
{
  if (numFRPConfinedConcrete == 0) {
    numFRPConfinedConcrete++;
    opserr
        << "FRPConfinedConcrete uniaxial material - Developed by Konstantinos "
           "G. Megalooikonomou University of Roma Tre Copyright 2009";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[18]; // size of arg list
  int numData;

  if (OPS_GetNumRemainingInputArgs() != 19) {
    opserr << "WARNING invalid #args: uniaxialMaterial FRPConfinedConcrete "
              "$tag $fpc1 $fpc2 $epsc0";
    opserr << " $D $c $Ej $Sj $tj $eju $S $fyl $fyh $dlong $dtrans $Es $v0 $k "
              "$useBuck\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial FRPConfinedConcrete tag"
           << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Material Properties: fpc1: Concrete Core "
              "Compressive Strength \n";
    opserr << "fpc2: Concrete Cover Compressive Strength \n";
    opserr << "epsc0: Strain Corresponding to Unconfined Concrete Strength \n";
    opserr << "D = Diameter of the Circular Section \n";
    opserr << "c = concrete cover \n";
    opserr << "Ej = Elastic Modulus of the Jacket \n";
    opserr
        << "Sj = Clear Spacing of the FRP strips - zero if it's continuous \n";
    opserr << "tj = Thickness of the FRP Jacket\n";
    opserr << "eju = Rupture strain of the Jacket\n";
    opserr << "S = Spacing of the stirrups\n";
    opserr << "fyl = Yielding Strength of longitudinal steel bars\n";
    opserr << "fyh = Yielding Strength of the hoops\n";
    opserr << "dlong = Diameter of the longitudinal bars\n";
    opserr << "dtrans = diameter of the stirrups\n";
    opserr << "Es = Steel's Elastic modulus\n";
    opserr << "vo = Poisson's coefficient for concrete\n";
    opserr
        << "k = reduction factor (0.5-0.8) for the rupture strain of the FRP\n";
    opserr << "useBuck = FRP Jacket Failure Criterion due to Buckling of "
              "Longitudinal Compressive Steel Bars (0 = not include it, 1= to "
              "include it)\n";
    return 0;
  }

  theMaterial = new FRPConfinedConcrete(
      iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
      dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
      dData[13], dData[14], dData[15], dData[16], dData[17]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "FRPConfinedConcrete\n";
    return 0;
  }

  return theMaterial;
}
