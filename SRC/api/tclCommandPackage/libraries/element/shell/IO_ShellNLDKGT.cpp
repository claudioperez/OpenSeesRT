

#include <g3_api.h>


#include <SRC/element/shell/ShellNLDKGT.h>
void *OPS_ShellNLDKGT(void)
{
  if (numShellNLDKGT == 0) {
    //    opserr << "Using ShellNLDKGT - Developed by:Shuhao Zhang & Xinzheng
    //    Lu";
    numShellNLDKGT++;
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellNLDKGT $tag $iNode $jNoe $kNode $secTag";
    return 0;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellNLDKGT \n";
    return 0;
  }

  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[4]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellNLDKGT " << iData[0] << "section "
           << iData[4] << " not found\n";
    return 0;
  }

  theElement =
      new ShellNLDKGT(iData[0], iData[1], iData[2], iData[3], *theSection);

  return theElement;
}
