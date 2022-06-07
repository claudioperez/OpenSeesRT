

#include <g3_api.h>


#include <SRC/element/UWelements/EmbeddedEPBeamInterface.h>
void *OPS_EmbeddedEPBeamInterface(void)
{
  if (num_EmbeddedEPBeamInterface == 0) {
    num_EmbeddedEPBeamInterface++;
    opserr << "EmbeddedEPBeamInterface element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedEPBeamInterface tag? \n";
    return 0;
  }

  int iData[1];
  int eleTag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedEPBeamInterface"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  theElement = new EmbeddedEPBeamInterface(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedEPBeamInterface\n";
    return 0;
  }

  return theElement;
}
