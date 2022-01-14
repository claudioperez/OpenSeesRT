

#include <g3_api.h>


#include <SRC/element/frictionBearing/FPBearingPTV.h>
void *OPS_FPBearingPTV()
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyBearing == 0) {
    opserr << "FPBearingPTV element - Written by Manish Kumar, University at "
              "Buffalo Copyright 2013\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new FPBearingPTV();
    return theEle;
  }

  if (numRemainingArgs < 30) {
    opserr << "ERROR - FPBearingPTV incorrect # args provided";
    return theEle;
  }

  int iData[13];
  double dData[17];
  int numData;

  // get the id and end nodes
  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];
  int iNode = iData[1];
  int jNode = iData[2];

  // MuRef
  double aDouble;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[0] = aDouble;
  double MuRef = dData[0];

  // IsPressureDependent
  numData = 1;
  int aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[3] = aInt;
  int IsPressureDependent = iData[3];

  // Reference Pressure
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[1] = aDouble;
  double pRef = dData[1];

  // IsTemperatureDependent
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[4] = aInt;
  int IsTemperatureDependent = iData[4];

  // Thermal diffusivity of Steel
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[2] = aDouble;
  double Diffusivity = dData[2];

  // Thermal conductivity of Steel
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[3] = aDouble;
  double Conductivity = dData[3];

  // IsVelocityDependent
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[5] = aInt;
  int IsVelocityDependent = iData[5];

  // rateParameter
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[4] = aDouble;
  double rateParameter = dData[4];

  // Effective radius of curvature of sliding surface
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[5] = aDouble;
  double ReffectiveFP = dData[5];

  // Radius of contact area
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[6] = aDouble;
  double Radius_Contact = dData[6];

  // kInitial
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[7] = aDouble;
  double kInitial = dData[7];

  // Material tags
  numData = 1;
  int aNumber;
  for (int i = 0; i < 4; i++) {
    if (OPS_GetIntInput(&numData, &aNumber) != 0) {
      opserr << "WARNING invalid material information\n";
      return 0;
    } else {
      iData[6 + i] = aNumber;
    }
  }

  int matIDP = iData[6];
  int matIDT = iData[7];
  int matIDMy = iData[8];
  int matIDMz = iData[9];

  UniaxialMaterial *theMaterialA = OPS_GetUniaxialMaterial(matIDP);
  UniaxialMaterial *theMaterialB = OPS_GetUniaxialMaterial(matIDT);
  UniaxialMaterial *theMaterialC = OPS_GetUniaxialMaterial(matIDMy);
  UniaxialMaterial *theMaterialD = OPS_GetUniaxialMaterial(matIDMz);

  // Orientation vector
  Vector x(3);
  Vector y(3);
  // x(0)=0.0; x(1)=0.0; x(2)=1.0;
  // y(0)=1.0; y(1)=0.0; y(2)=0.0;
  numData = 1;
  for (int i = 0; i < 6; i++) {
    if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
      opserr << "WARNING invalid element data\n";
      return 0;
    } else {
      dData[8 + i] = aDouble;
    }
  }
  x(0) = dData[8];
  x(1) = dData[9];
  x(2) = dData[10];
  y(0) = dData[11];
  y(1) = dData[12];
  y(2) = dData[13];

  // Shear distance
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[14] = aDouble;
  double shearDist = dData[14];

  // doRayleigh
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[10] = aInt;
  int doRayleigh = iData[10];

  // mass
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[15] = aDouble;
  double mass = dData[15];

  // NumberOfIterations
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[11] = aInt;
  int iter = iData[11];

  // tol
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag
           << endln;
    return 0;
  }
  dData[16] = aDouble;
  double tol = dData[16];

  // Units 1: N,m,s,C; 2: kN,m,s,C; 3: N,mm,s,C; 4: kN,mm,s,C; 5: lb,in,s,C; 6:
  // kip,in,s,C; 7: lb,ft,s,C; 8: kip,ft,s,C
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[12] = aInt;
  int unit = iData[12];

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm == 3) {

    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - FPBearingPTV \n";
    }

    theEle = new FPBearingPTV(
        eleTag, iNode, jNode, MuRef, IsPressureDependent, pRef,
        IsTemperatureDependent, Diffusivity, Conductivity, IsVelocityDependent,
        rateParameter, ReffectiveFP, Radius_Contact, kInitial, *theMaterialA,
        *theMaterialB, *theMaterialC, *theMaterialD, x, y, shearDist,
        doRayleigh, mass, iter, tol, unit);
  }

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag
           << endln;
    delete theMaterialA;
    delete theMaterialB;
    delete theMaterialC;
    delete theMaterialD;
    return 0;
  }

  return theEle;
}
