

#include <g3_api.h>


#include <SRC/material/uniaxial/UVCuniaxial.h>
void *OPS_UVCuniaxial(void)
{
  if (numUVCuniaxial == 0) {
    opserr << "Using the UVCuniaxial material, see "
              "https://www.epfl.ch/labs/resslab/resslab-tools/"
           << endln;
    numUVCuniaxial++;
  }
  UniaxialMaterial *theMaterial = 0;

  // Parameters for parsing
  const int N_TAGS = 1;
  const int N_BASIC_PROPERTIES = 4;
  const int N_UPDATED_PROPERTIES = 2;
  const int N_PARAM_PER_BACK = 2;
  const int MAX_BACKSTRESSES = 8;
  const int BACKSTRESS_SPACE = MAX_BACKSTRESSES * N_PARAM_PER_BACK;

  std::string inputInstructions =
      "Invalid args, want:\n"
      "uniaxialMaterial UVCuniaxial "
      "tag? E? fy? QInf? b? DInf? a? "
      "N? C1? gamma1? <C2? gamma2? C3? gamma3? ... C8? gamma8?>\n"
      "Note: to neglect the updated model, set DInf = 0.0";

  // Containers for the inputs
  int nInputsToRead;
  int nBackstresses[1];                     // for N
  int materialTag[1];                       // for the tag
  double basicProps[4];                     // holds E, fy, QInf, b
  double updProps[2];                       // holds DInf, a
  double backstressProps[BACKSTRESS_SPACE]; // holds C's and gamma's
  std::vector<double> cK;
  std::vector<double> gammaK;

  // Get the material tag
  nInputsToRead = N_TAGS;
  if (OPS_GetIntInput(&nInputsToRead, materialTag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial UVCuniaxial tag" << endln;
    return 0;
  }

  // Get E, fy, qInf, b
  nInputsToRead = N_BASIC_PROPERTIES;
  if (OPS_GetDoubleInput(&nInputsToRead, basicProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }

  // Read in the updated model paramters
  nInputsToRead = N_UPDATED_PROPERTIES;
  if (OPS_GetDoubleInput(&nInputsToRead, updProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }

  // Get the number of backstresses
  nInputsToRead = 1;
  if (OPS_GetIntInput(&nInputsToRead, nBackstresses) != 0) {
    opserr << "WARNING N must be an integer" << inputInstructions.c_str()
           << endln;
    return 0;
  }
  if (nBackstresses[0] > MAX_BACKSTRESSES) {
    opserr << "WARNING: Too many backstresses defined, maximum is: "
           << MAX_BACKSTRESSES << endln << inputInstructions.c_str() << endln;
    return 0;
  }

  // Get the backstress parameters
  nInputsToRead = 2 * nBackstresses[0];
  if (OPS_GetDoubleInput(&nInputsToRead, backstressProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }
  // cK's alternate with gammaK's
  for (int i = 0; i < nBackstresses[0]; ++i) {
    cK.push_back(backstressProps[2 * i]);
    gammaK.push_back(backstressProps[1 + 2 * i]);
  }

  // Allocate the material
  theMaterial = new UVCuniaxial(materialTag[0], basicProps[0], basicProps[1],
                                basicProps[2], basicProps[3], updProps[0],
                                updProps[1], cK, gammaK);

  return theMaterial;
}
