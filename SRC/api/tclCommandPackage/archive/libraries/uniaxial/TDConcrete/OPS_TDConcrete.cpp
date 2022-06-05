

#include <g3_api.h>


#include <SRC/material/uniaxial/TDConcrete.h>
void *OPS_TDConcrete(G3_Runtime *rt)
{
  // Print description of material model:
  if (numTDConcrete == 0) {
    opserr << "Time-Dependent Concrete Material Model - Written by Adam "
              "Knaack, University of Notre Dame, 2012 \n";
    numTDConcrete = 1;
  }

  // Pointer to a uniaxial material that will be returned:
  UniaxialMaterial *theMaterial = 0;

  // Parse the input line for the material parameters:
  int iData;
  int numData;
  int numArgs;

  numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs == 13) {
    // TDConcrete(int tag, double _fc, double _epsc0, double _fcu,
    // double _epscu, double _tcr, double _ft, double _Ets, double _Ec, double
    // _age, double _epsshu)
    double dData[12];

    // Collect material tag:
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial TDConcrete tag\n";
      return 0;
    }

    // Collect input parameters:
    numData = 12;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING: invalid material property definition\n";
      return 0;
    }

    // Create a new materiadouble
    theMaterial = new TDConcrete(iData, dData[0], dData[1], dData[2], dData[3],
                                 dData[4], dData[5], dData[6], dData[7],
                                 dData[8], dData[9], dData[10], dData[11]);
    if (theMaterial == 0) {
      opserr
          << "WARNING: could not create uniaxialMaterial of type TDConcrete \n";
      return 0;
    }

    // Return new material:
    return theMaterial;
  }
}
