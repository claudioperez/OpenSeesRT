

#include <g3_api.h>


#include <SRC/element/shell/ShellDKGQ.h>
void *OPS_ShellDKGQ(void)
{
  if (numShellDKGQ == 0) {
    //    opserr << "Using ShellDKGQ - Developed by: Lisha Wang and Xinzheng
    //    Lu\n";
    numShellDKGQ++;
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKGQ \n";
    return 0;
  }

  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellDKGQ " << iData[0] << "section " << iData[5]
           << " not found\n";
    return 0;
  }

  theElement = new ShellDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                             *theSection);

  return theElement;
}
