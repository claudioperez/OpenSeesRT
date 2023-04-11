#include <assert.h>
#include <tcl.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ShellDKGQ.h>

Element*
TclDispatch_newShellDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  // if (numShellDKGQ == 0) {
  //   //    opserr << "Using ShellDKGQ - Developed by: Lisha Wang and Xinzheng
  //   //    Lu\n";
  //   numShellDKGQ++;
  // }

  Element *theElement = 0;

  if (argc < 6) {
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
      builder->getSection(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellDKGQ " << iData[0] << "section " << iData[5]
           << " not found\n";
    return 0;
  }

  theElement = new ShellDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                             *theSection);

  return theElement;
}
