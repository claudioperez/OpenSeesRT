

#include <g3_api.h>


#include <SRC/element/shell/ASDShellQ4.h>
void *OPS_ASDShellQ4(void)
{
  static bool first_done = false;
  if (!first_done) {
    opserr << "Using ASDShellQ4 - Developed by: Massimo Petracca, Guido "
              "Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 6) {
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

  if (numArgs == 7) {
    const char *type = OPS_GetString();
    if ((strcmp(type, "-corotational") == 0) ||
        (strcmp(type, "-Corotational") == 0))
      corotational = true;
  }

  SectionForceDeformation *section = OPS_getSectionForceDeformation(iData[5]);

  if (section == 0) {
    opserr << "ERROR:  element ASDShellQ4 " << iData[0] << "section "
           << iData[5] << " not found\n";
    return 0;
  }

  return new ASDShellQ4(iData[0], iData[1], iData[2], iData[3], iData[4],
                        section, corotational);
}
