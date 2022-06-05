

#include <g3_api.h>


#include <SRC/material/uniaxial/IMKPeakOriented.h>
void *OPS_IMKPeakOriented(G3_Runtime *rt)
{
  if (numIMKPeakOrientedMaterials == 0) {
    numIMKPeakOrientedMaterials++;
    OPS_Error("IMK Model with Peak-Oriented Response - Code by A. ELKADY & H. "
              "ELJISR (July 2020)\n",
              1);
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[23];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial IMKPeakOriented tag" << endln;
    return 0;
  }

  numData = 23;

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial IMKPeakOriented tag? Ke? ";
    opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FmaxFy_pos? ResF_pos? ";
    opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FmaxFy_neg? ResF_neg? ";
    opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new IMKPeakOriented(
      iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
      dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
      dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
      dData[19], dData[20], dData[21], dData[22]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "IMKPeakOriented Material\n";
    return 0;
  }

  return theMaterial;
}
