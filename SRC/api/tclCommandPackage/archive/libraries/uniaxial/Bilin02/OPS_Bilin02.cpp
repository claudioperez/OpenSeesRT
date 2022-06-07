

#include <g3_api.h>


#include <SRC/material/uniaxial/Bilin02.h>
void *OPS_Bilin02(void)
{
  if (numBilin02Materials == 0) {
    numBilin02Materials++;
    opserr
        << "Modified Ibarra-Medina-Krawinkler Model with Bilinear Hysteretic "
           "Response\n"; // Updated: Filipe Ribeiro and Andre Barbosa
    opserr
        << "Implementation and Calibration for CPH and FLPH by F.L.A. Ribeiro "
           "and A.R. Barbosa\n"; // Updated: Filipe Ribeiro and Andre Barbosa
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[24]; // Updated: Filipe Ribeiro and Andre Barbosa
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Bilin02 tag" << endln;
    return 0;
  }

  // Changed in order to account for the nFactor as an optional input    //
  // Updated: Filipe Ribeiro and Andre Barbosa numData = 24;               //
  // Updated: Filipe Ribeiro and Andre Barbosa
  numData = OPS_GetNumRemainingInputArgs(); // Updated: Filipe Ribeiro and Andre
                                            // Barbosa

  if (numData != 23 &&
      numData != 24) { // Updated: Filipe Ribeiro and Andre Barbosa
    opserr << "Invalid Args want: uniaxialMaterial Bilin02 tag? Ke? AsPos? "
              "AsNeg? My_pos? My_neg? LamdaS? ";
    opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? "
              "Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? ";
    opserr
        << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg?  <nFactor?> \n"; // Updated:
                                                                           // Filipe
                                                                           // Ribeiro
                                                                           // and
                                                                           // Andre
                                                                           // Barbosa
    return 0; // Updated: Filipe Ribeiro and Andre Barbosa
  }

  //  if (OPS_GetDoubleInput(&numData, dData) != 0) {                   //
  //  Updated: Filipe Ribeiro and Andre Barbosa
  //    opserr << "Invalid Args want: uniaxialMaterial Bilin02 tag? Ke? nFactor?
  //    AsPos? AsNeg? My_pos? My_neg? LamdaS? "; opserr << "LamdaD?  LamdaA?
  //    LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? Thetap_neg? Thetapc_pos?
  //    Thetapc_neg?KPos? "; opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?
  //    PDNeg\n"; return 0;
  //  } // Updated: Filipe Ribeiro and Andre Barbosa

  if (numData == 23) {
    if (OPS_GetDoubleInput(&numData, dData) !=
        0) { // Updated: Filipe Ribeiro and Andre Barbosa
      opserr << "Invalid Args want: uniaxialMaterial Bilin02 tag? Ke? AsPos? "
                "AsNeg? My_pos? My_neg? LamdaS? ";
      opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? "
                "Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? ";
      opserr
          << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg? <nFactor?> \n"; // Updated: Filipe Ribeiro and Andre Barbosa
      return 0;
    }
    // Parsing was successful, allocate the material
    theMaterial = new Bilin02(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
        dData[19], dData[20], dData[21],
        dData[22]); // Updated: Filipe Ribeiro and Andre Barbosa

  } else if (numData == 24) { // Updated: Filipe Ribeiro and Andre Barbosa
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial Bilin02 tag? Ke? AsPos? "
                "AsNeg? My_pos? My_neg? LamdaS? ";
      opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? "
                "Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? ";
      opserr
          << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg? <nFactor?>\n"; // Updated:
                                                                           // Filipe
                                                                           // Ribeiro
                                                                           // and
                                                                           // Andre
                                                                           // Barbosa
      return 0;
    }
    // Parsing was successful, allocate the material
    theMaterial = new Bilin02(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18],
        dData[19], dData[20], dData[21], dData[22],
        dData[23]); // Updated: Filipe Ribeiro and Andre Barbosa
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Bilin02 "
              "Material\n";
    return 0;
  }

  return theMaterial;
}
