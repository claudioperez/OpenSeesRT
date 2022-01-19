

#include <g3_api.h>


#include <SRC/element/shell/ShellANDeS.h>
void *OPS_ShellANDeS(void)
{

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellANDeS $tag $iNode $jNode $kNode $thick $E "
              "$nu $rho";
    return 0;
  }

  int iData[4];
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellANDeS \n";
    return 0;
  }

  double dData[11];
  numArgs = OPS_GetNumRemainingInputArgs();
  if (OPS_GetDoubleInput(&numArgs, dData) != 0) {
    opserr << "WARNING invalid double thickness: element ShellANDeS \n";
    return 0;
  }
  if (numArgs == 4) {
    theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3],
                                dData[0], dData[1], dData[2], dData[3]);
  } else if (numArgs == 11) {
    theElement =
        new ShellANDeS(iData[0], iData[1], iData[2], iData[3], dData[0],
                       dData[1], dData[2], dData[3], dData[4], dData[5],
                       dData[6], dData[7], dData[8], dData[9], dData[10]);
  }

  return theElement;
}
