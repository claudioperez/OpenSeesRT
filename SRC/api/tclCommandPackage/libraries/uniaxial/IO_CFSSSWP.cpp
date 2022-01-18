

#include <g3_api.h>


#include <SRC/material/uniaxial/CFSSSWP.h>
void *OPS_CFSSSWP(void)
{
  // print out some KUDO's
  if (numCFSSSWP == 0) {
    opserr << "Cold Formed Steel Steel-Sheathed Shear Wall Panel "
              "uniaxialMaterial - Written by Smail KECHIDI Ph.D Student at "
              "University of Blida 1 - Please when using this make reference "
              "as: Smail Kechidi and Nouredine Bourahla (2016), Deteriorating "
              "hysteresis model for cold-formed steel shear wall panel based "
              "on its physical and mechanical characteristics, Journal of "
              "Thin-Walled Structures, DOI: 10.1016/j.tws.2015.09.022\n";
    numCFSSSWP = 1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int iData[1];
  double dData[16];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial CFSSSWP tag" << endln;
    return 0;
  }

  numData = 15;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Material parameters\n";
    return 0;
  }

  //
  // create a new material
  //

  theMaterial =
      new CFSSSWP(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                  dData[5], dData[6], dData[7], dData[8], dData[9], dData[10],
                  dData[11], dData[12], dData[13], dData[14]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type CFSSSWP\n";
    return 0;
  }

  // return the material
  return theMaterial;
}
