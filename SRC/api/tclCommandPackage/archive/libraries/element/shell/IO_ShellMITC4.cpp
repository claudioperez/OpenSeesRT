

#include <g3_api.h>


#include <SRC/element/shell/ShellMITC4.h>
void *OPS_ShellMITC4(void)
{
  if (numShellMITC4 == 0) {
    //    opserr << "Using ShellMITC4 - Developed by: Leopoldo Tesser, Diego A.
    //    Talledo, Véronique Le Corvec\n";
    numShellMITC4++;
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4 $tag $iNode $jNoe $kNode $lNode "
              "$secTag<-updateBasis>";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC4 \n";
    return 0;
  }
  bool updateBasis = false;

  if (numArgs == 7) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-updateBasis") == 0)
      updateBasis = true;
  }

  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC4 " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  theElement = new ShellMITC4(iData[0], iData[1], iData[2], iData[3], iData[4],
                              *theSection, updateBasis);

  return theElement;
}
