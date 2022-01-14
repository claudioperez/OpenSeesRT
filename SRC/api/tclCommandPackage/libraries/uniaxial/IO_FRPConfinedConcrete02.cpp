

#include <g3_api.h>


#include <SRC/material/uniaxial/FRPConfinedConcrete02.h>
void *OPS_FRPConfinedConcrete02(G3_Runtime *rt)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial FRPConfinedConcrete02 tag"
           << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs(); // get number of input args
                                            // excluding the material tag

  if (numData != 6 && numData != 9 && numData != 11) {
    opserr << "Incorrect # args, want: uniaxialMaterial FRPConfinedConcrete02 "
              "tag? fc0? ec0? Ec? ft? Ets? Unit?"
           << endln;
    opserr << "Or: uniaxialMaterial FRPConfinedConcrete02 tag? fc0? ec0? Ec? "
              "-Ultimate fcc? ecu? ft? Ets? Unit?"
           << endln;
    opserr << "Or: uniaxialMaterial FRPConfinedConcrete02 tag? fc0? ec0? Ec? "
              "-JacketC t? Efrp? eps_h_rup? R? ft? Ets? Unit?"
           << endln;
    return 0;
  }

  if (numData == 6) { // Input for unconfined concrete
    int numData1 = 6;
    double dData[6];
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 "
             << iData[0] << "fc0? ec0? Ec? ft? Ets? Unit?" << endln;
      return 0;
    }
    theMaterial =
        new FRPConfinedConcrete02(iData[0], dData[0], dData[1], dData[2],
                                  dData[3], dData[4], (int)dData[5]);

  } else if (numData == 9) { // Ultimate stress/strain input by users
    double dData[8];
    int numData1 = 3;
    int numData2 = 5;
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 "
             << iData[0] << "fc0? ec0? Ec? -Ultimate fcc? ecu? ft? Ets? Unit?"
             << endln;
      return 0;
    }

    const char *str = OPS_GetString();
    // OPS_GetStringCopy(&str);
    if (strcmp(str, "-Ultimate") == 0) {
      if (OPS_GetDoubleInput(&numData2, dData + 3) != 0) {
        opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 "
               << iData[0] << "fc0? ec0? Ec? -Ultimate fcc? ecu? ft? Ets? Unit?"
               << endln;
        return 0;
      }
    } else {
      opserr << "Invalid input parameter for uniaxialMaterial "
                "FRPConfinedConcrete02 with tag  "
             << iData[0] << ", want: -Ultimate" << endln;
      return 0;
    }
    theMaterial = new FRPConfinedConcrete02(iData[0], dData[0], dData[1],
                                            dData[2], dData[3], dData[4],
                                            dData[5], dData[6], (int)dData[7]);
  } else { // FRP-confined concrete in circular columns, FRP jacket properties
           // input by users
    double dData[10];
    int numData1 = 3;
    int numData2 = 7;
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 "
             << iData[0]
             << "fc0? ec0? Ec? -JacketC tfrp? Efrp? erup? R? ft? Ets? Unit?"
             << endln;
      return 0;
    }

    const char *str = OPS_GetString();
    // OPS_GetStringCopy(&str);
    if (strcmp(str, "-JacketC") == 0) {
      if (OPS_GetDoubleInput(&numData2, dData + 3) != 0) {
        opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 "
               << iData[0]
               << "fc0? ec0? Ec? -JacketC tfrp? Efrp? erup? R? ft? Ets? Unit?"
               << endln;
        return 0;
      }
    } else {
      opserr << "Invalid input parameter for uniaxialMaterial "
                "FRPConfinedConcrete02 with tag "
             << iData[0] << ", want: -JacketC" << endln;
      return 0;
    }
    theMaterial = new FRPConfinedConcrete02(
        iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
        dData[6], dData[7], dData[8], (int)dData[9]);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial FRPConfinedConcrete02 "
           << iData[0] << endln;
    return 0;
  }

  return theMaterial;
}
