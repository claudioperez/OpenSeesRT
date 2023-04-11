#include <assert.h>
#include <tcl.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ASDShellQ4.h>

Element*
TclDispatch_newASDShellQ4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  // static bool first_done = false;
  // if (!first_done) {
  //   opserr << "Using ASDShellQ4 - Developed by: Massimo Petracca, Guido "
  //             "Camata, ASDEA Software Technology\n";
  //   first_done = true;
  // }

  if (argc < 6) {
    opserr << "Want: element ASDShellQ4 $tag $iNode $jNode $kNode $lNode "
              "$secTag <-corotational>";
    return 0;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ASDShellQ4 \n";
    return 0;
  }
  bool corotational = false;

  if (argc == 7) {
    const char *type = OPS_GetString();
    if ((strcmp(type, "-corotational") == 0) ||
        (strcmp(type, "-Corotational") == 0))
      corotational = true;
  }

  SectionForceDeformation *section = builder->getSection(iData[5]);

  if (section == 0) {
    opserr << "ERROR:  element ASDShellQ4 " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  return new ASDShellQ4(iData[0], iData[1], iData[2], iData[3], iData[4],
                        section, corotational);
}
