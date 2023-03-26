
#include <assert.h>

#include <tcl.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>


#include <element/Shell/ShellNLDKGQ.h>

Element*
TclDispatch_newShellNLDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // if (numShellNLDKGQ == 0) {
  //   //    opserr << "Using ShellNLDKGQ - Developed by: Lisha Wang,Xinzheng Lu
  //   //    and Quan Gu\n";
  //   numShellNLDKGQ++;
  // }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr
        << "Want: element ShellNLDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellNLDKGQ \n";
    return 0;
  }

  SectionForceDeformation *theSection =
      builder->getSection(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellNLDKGQ " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  theElement = new ShellNLDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                               *theSection);

  return theElement;
}
