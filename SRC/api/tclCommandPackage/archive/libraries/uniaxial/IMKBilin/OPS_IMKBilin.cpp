

#include <g3_api.h>


#include <SRC/material/uniaxial/IMKBilin.h>
void *OPS_IMKBilin(void)
{
  if (numIMKBilinMaterials == 0) {
    numIMKBilinMaterials++;
    OPS_Error("Mod. IMK Model with Bilinear Hysteretic Response - Code by "
              "A.ELKADY-21\n",
              1);
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[21];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial IMKBilin tag" << endln;
    return 0;
  }

  numData = 21;

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial IMKBilin tag? Ke? ";
    opserr << "Theta_p_pos? Theta_pc_pos? Theta_u_pos? Mpe_pos? MmaxMpe_pos? "
              "ResM_pos? ";
    opserr << "Theta_p_neg? Theta_pc_neg? Theta_u_neg? Mpe_neg? MmaxMpe_neg? "
              "ResM_neg? ";
    opserr << "LamdaS?  LamdaC? LamdaK? Cs? Cc? Ck? D_pos? D_neg? ";
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new IMKBilin(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
                   dData[5], dData[6], dData[7], dData[8], dData[9], dData[10],
                   dData[11], dData[12], dData[13], dData[14], dData[15],
                   dData[16], dData[17], dData[18], dData[19], dData[20]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type IMKBilin "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
