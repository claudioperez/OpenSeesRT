

#include <g3_api.h>


#include <SRC/element/UWelements/EmbeddedBeamInterfaceL.h>
void *OPS_EmbeddedBeamInterfaceL(void)
{
  if (num_EmbeddedBeamInterfaceL == 0) {
    num_EmbeddedBeamInterfaceL++;
    opserr << "EmbeddedBeamInterfaceL element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedBeamInterfaceL tag? \n";
    return 0;
  }

  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceL"
           << endln;
    return 0;
  }

  int eleTag = iData[0];

  theElement = new EmbeddedBeamInterfaceL(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedBeamInterfaceL\n";
    return 0;
  }

  return theElement;
}
