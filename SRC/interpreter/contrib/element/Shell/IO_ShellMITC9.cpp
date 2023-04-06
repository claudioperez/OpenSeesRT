
#include <assert.h>

#include <tcl.h>
#include <BasicModelBuilder.h>


#include <element/Shell/ShellMITC9.h>

Element*
TclDispatch_newShellMITC9(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // if (numShellMITC9 == 0) {
  //   opserr << "Using ShellMITC9 - Developed by: Leopoldo Tesser and Diego A. "
  //             "Talledo\n";
  //   numShellMITC9++;
  // }

  Element *theElement = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 11) {
    opserr << "Want: element ShellMITC9 $tag $node1 $node2 .... $node9 $secTag";
    return 0;
  }

  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC9\n";
    return 0;
  }

  SectionForceDeformation *theSection =
      builder->getSection(iData[10]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC9 " << iData[0] << "section "
           << iData[10] << " not found\n";
    return 0;
  }

  theElement =
      new ShellMITC9(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5],
                     iData[6], iData[7], iData[8], iData[9], *theSection);

  return theElement;
}
