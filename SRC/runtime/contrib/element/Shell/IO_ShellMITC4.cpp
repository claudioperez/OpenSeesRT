
#include <assert.h>

#include <tcl.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ShellMITC4.h>

Element*
TclDispatch_newShellMITC4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // if (numShellMITC4 == 0) {
  //   //    opserr << "Using ShellMITC4 - Developed by: Leopoldo Tesser, Diego A.
  //   //    Talledo, V�ronique Le Corvec\n";
  //   numShellMITC4++;
  // }

  Element *theElement = nullptr;

  if (argc < 6) {
    opserr << "Want: element ShellMITC4 $tag $iNode $jNode $kNode $lNode "
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

  if (argc == 7) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-updateBasis") == 0)
      updateBasis = true;
  }

  SectionForceDeformation *theSection =
      builder->getSection(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC4 " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  theElement = new ShellMITC4(iData[0], iData[1], iData[2], iData[3], iData[4],
                              *theSection, updateBasis);

  return theElement;
}
