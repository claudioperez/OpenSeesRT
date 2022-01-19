

#include <g3_api.h>


#include <SRC/material/uniaxial/IMKPinching.h>
void *OPS_IMKPinching(G3_Runtime *rt)
{
  if (numIMKPinchingMaterials == 0) {
    numIMKPinchingMaterials++;
    OPS_Error("IMK Model with Pinched Response - Code by A. ELKADY & H. ELJISR "
              "(July 2020)\n",
              1);
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[25];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial IMKPinching tag" << endln;
    return 0;
  }

  numData = 25;

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial IMKPinching tag? Ke? ";
    opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FmaxFy_pos? ResF_pos? ";
    opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FmaxFy_neg? ResF_neg? ";
    opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? "
              "kappaF? kappaD? ";
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new IMKPinching(
      iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
      dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
      dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
      dData[19], dData[20], dData[21], dData[22], dData[23], dData[24]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type IMKPinching "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
