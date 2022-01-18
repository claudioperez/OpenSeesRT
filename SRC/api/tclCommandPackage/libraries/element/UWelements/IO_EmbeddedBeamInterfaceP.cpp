

#include <g3_api.h>


#include <SRC/element/UWelements/EmbeddedBeamInterfaceP.h>
void *OPS_EmbeddedBeamInterfaceP(void)
{
  if (num_EmbeddedBeamInterfaceP == 0) {
    num_EmbeddedBeamInterfaceP++;
    opserr << "EmbeddedBeamInterfaceP element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedBeamInterfaceP tag? \n";
    return 0;
  }

  int iData[1];
  int eleTag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceP"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  theElement = new EmbeddedBeamInterfaceP(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedBeamInterfaceP\n";
    return 0;
  }

  return theElement;
}
