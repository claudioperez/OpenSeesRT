
#include <assert.h>

#include <tcl.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ShellMITC4Thermal.h>

Element*
TclDispatch_newShellMITC4Thermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // if (numShellMITC4Thermal == 0) {
  //   opserr << "Using ShellMITC4Thermal - Developed by: Leopoldo Tesser, Diego "
  //             "A. Talledo, Veronique Le Corvec\n";
  //   numShellMITC4Thermal++;
  // }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4Thermal $tag $iNode $jNoe $kNode $lNode "
              "$secTag";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC4Thermal \n";
    return 0;
  }

  SectionForceDeformation *theSection =
      builder->getSection(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC4Thermal " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  theElement = new ShellMITC4Thermal(iData[0], iData[1], iData[2], iData[3],
                                     iData[4], *theSection);

  return theElement;
}
