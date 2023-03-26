
#include <assert.h>

#include <tcl.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ShellNLDKGT.h>

Element*
TclDispatch_newShellNLDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // if (numShellNLDKGT == 0) {
  //   //    opserr << "Using ShellNLDKGT - Developed by:Shuhao Zhang & Xinzheng
  //   //    Lu";
  //   numShellNLDKGT++;
  // }

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
      builder->getSection(iData[4]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellNLDKGT " << iData[0] << "section "
           << iData[4] << " not found\n";
    return 0;
  }

  theElement =
      new ShellNLDKGT(iData[0], iData[1], iData[2], iData[3], *theSection);

  return theElement;
}
