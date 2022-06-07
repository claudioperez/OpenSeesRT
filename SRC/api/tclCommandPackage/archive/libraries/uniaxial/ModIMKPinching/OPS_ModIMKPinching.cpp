

#include <g3_api.h>


#include <SRC/material/uniaxial/ModIMKPinching.h>
void *OPS_ModIMKPinching(G3_Runtime *rt)
{
  if (numModIMKPinchingMaterials == 0) {
    numModIMKPinchingMaterials++;
    opserr << "Modified Ibarra-Medina-Krawinkler Model with Pinched Hysteretic "
              "Response\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[27]; // Updated: Filipe Ribeiro and Andre Barbosa
  int numData = 1;
  // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ModIMKPinching tag" << endln;
    return 0;
  }

  // Changed in order to account for the nFactor as an optional input    //
  // Updated: Filipe Ribeiro and Andre Barbosa numData = 27;

  numData = OPS_GetNumRemainingInputArgs(); // Updated: Filipe Ribeiro and Andre
                                            // Barbosa
  if (numData != 27 && numData != 26) {
    opserr << "Invalid Args want: uniaxialMaterial ModIMKPinching tag? Ke?, "
              "alfaPos?, alfaNeg?, My_pos?, My_neg?";
    opserr << "FprPos?, FprNeg?, A_pinch?, Ls?, Ld?, La?, Lk?, Cs?, Cd?, Ca?, "
              "Ck?, thetaPpos?, thetaPneg?";
    opserr
        << "thetaPCpos?, thetaPCneg?, ResfacPos?, ResfacNeg?, fracDispPos?, "
           "fracDispNeg?,DPos?, DNeg?, <nFactor?>"; // Updated: Filipe Ribeiro
                                                    // and Andre Barbosa

    return 0;
  }

  if (numData == 26) { // Updated: Filipe Ribeiro and Andre Barbosa
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial ModIMKPinching tag? Ke?, "
                "alfaPos?, alfaNeg?, My_pos?, My_neg?";
      opserr << "FprPos?, FprNeg?, A_pinch?, Ls?, Ld?, La?, Lk?, Cs?, Cd?, "
                "Ca?, Ck?, thetaPpos?, thetaPneg?";
      opserr
          << "thetaPCpos?, thetaPCneg?, ResfacPos?, ResfacNeg?, fracDispPos?, "
             "fracDispNeg?,DPos?, DNeg?, <nFactor?>"; // Updated: Filipe Ribeiro
                                                      // and Andre Barbosa

      return 0;
    }

    // Parsing was successful, allocate the material with zero index
    theMaterial = new ModIMKPinching(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
        dData[19], dData[20], dData[21], dData[22], dData[23], dData[24],
        dData[25]);

  } else if (numData == 27) { // Updated: Filipe Ribeiro and Andre Barbosa
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial ModIMKPinching tag? Ke?, "
                "alfaPos?, alfaNeg?, My_pos?, My_neg?";
      opserr << "FprPos?, FprNeg?, A_pinch?, Ls?, Ld?, La?, Lk?, Cs?, Cd?, "
                "Ca?, Ck?, thetaPpos?, thetaPneg?";
      opserr
          << "thetaPCpos?, thetaPCneg?, ResfacPos?, ResfacNeg?, fracDispPos?, "
             "fracDispNeg?,DPos?, DNeg?, <nFactor?>"; // Updated: Filipe Ribeiro
                                                      // and Andre Barbosa

      return 0;
    }

    // Parsing was successful, allocate the material with zero index
    theMaterial = new ModIMKPinching(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
        dData[19], dData[20], dData[21], dData[22], dData[23], dData[24],
        dData[25], dData[26]); // Updated: Filipe Ribeiro and Andre Barbosa
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ModIMKPinching Material\n";
    return 0;
  }

  return theMaterial;
}
