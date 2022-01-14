

#include <g3_api.h>


#include <SRC/material/uniaxial/backbone/MultilinearBackbone.h>
void *OPS_MultilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Invalid number of args, want: hystereticBackbone Multilinear "
              "tag? e1? s1? e2? s2? ..."
           << endln;
    return 0;
  }

  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Multilinear" << endln;
    return 0;
  }

  int numPoints = OPS_GetNumRemainingInputArgs() / 2;
  numData = 2 * numPoints;
  Vector e(numPoints);
  Vector s(numPoints);
  double *dData = new double[numData];
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Multilinear"
           << endln;
    return 0;
  }
  for (int i = 0; i < numPoints; i++) {
    e(i) = dData[2 * i];
    s(i) = dData[2 * i + 1];
  }

  theBackbone = new MultilinearBackbone(iData[0], numPoints, e, s);
  if (theBackbone == 0) {
    opserr << "WARNING could not create MultilinearBackbone\n";
    return 0;
  }

  delete[] dData;

  return theBackbone;
}
