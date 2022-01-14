

#include <g3_api.h>


#include <SRC/material/uniaxial/backbone/ManderBackbone.h>
void *OPS_ManderBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "Invalid number of args, want: hystereticBackbone Mander tag? "
              "fc? epsc? E?"
           << endln;
    return 0;
  }

  int iData[1];
  double dData[3];

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Mander" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Mander" << endln;
    return 0;
  }

  theBackbone = new ManderBackbone(iData[0], dData[0], dData[1], dData[2]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create ManderBackbone\n";
    return 0;
  }

  return theBackbone;
}
