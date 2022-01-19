

#include <g3_api.h>


#include <SRC/element/UWelements/QuadBeamEmbedContact.h>
void *OPS_QuadBeamEmbedContact(void)
{
  if (num_QuadBeamEmbedContact == 0) {
    num_QuadBeamEmbedContact++;
    opserr << "QuadBeamEmbedContact element - Written: A.Ghofrani, P.Arduino, "
              "U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 10) {
    opserr << "Want: QuadBeamEmbedContact tag? Qnd1? Qnd2? Qnd3? Qnd4? Bnd1? "
              "Bnd2? radius? fricCoeff? normalPenalty? <tangentialPenalty?> \n";
    return 0;
  }

  int iData[7];
  int eleTag = 0;
  int numData = 7;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element QuadBeamEmbedContact"
           << endln;
    return 0;
  }

  numData = 3;
  double dData[3];
  double oData[1];

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
    return 0;
  }

  oData[0] = dData[2];
  numData = numArgs - 10;
  if (numData != 0)
    if (OPS_GetDouble(&numData, oData) != 0) {
      opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
      return 0;
    }

  eleTag = iData[0];

  theElement = new QuadBeamEmbedContact(iData[0], iData[1], iData[2], iData[3],
                                        iData[4], iData[5], iData[6], dData[0],
                                        dData[1], dData[2], oData[0]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type QuadBeamEmbedContact\n";
    return 0;
  }

  return theElement;
}
