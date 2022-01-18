

#include <g3_api.h>


#include <SRC/material/uniaxial/backbone/TrilinearBackbone.h>
void *OPS_TrilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Invalid number of args, want: hystereticBackbone Trilinear tag? "
              "e1? s1? e2? s2? e3? s3?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[6];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Trilinear" << endln;
    return 0;
  }

  numData = 6;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Trilinear" << endln;
    return 0;
  }

  theBackbone = new TrilinearBackbone(iData[0], dData[0], dData[1], dData[2],
                                      dData[3], dData[4], dData[5]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create TrilinearBackbone\n";
    return 0;
  }

  return theBackbone;
}
void *OPS_BilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Invalid number of args, want: hystereticBackbone Bilinear tag? "
              "e1? s1? e2? s2?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[4];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Bilinear" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Bilinear" << endln;
    return 0;
  }

  theBackbone =
      new TrilinearBackbone(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create BilinearBackbone\n";
    return 0;
  }

  return theBackbone;
}
