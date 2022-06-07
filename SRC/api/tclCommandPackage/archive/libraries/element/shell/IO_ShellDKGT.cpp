

#include <g3_api.h>


#include <SRC/element/shell/ShellDKGT.h>
void *OPS_ShellDKGT(void)
{
  if (numShellDKGT == 0) {
    //    opserr << "Using ShellDKGT - Developed by: Shuhao Zhang and Xinzheng
    //    Lu\n";
    numShellDKGT++;
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellDKGT $tag $iNode $jNoe $kNode $secTag";
    return 0;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKGT \n";
    return 0;
  }

  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[4]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellDKGT " << iData[0] << "section " << iData[4]
           << " not found\n";
    return 0;
  }

  double b_data[3] = {0, 0, 0};

  int num_remaining_args = OPS_GetNumRemainingInputArgs();
  if (num_remaining_args > 3) {
    num_remaining_args = 3;
  }
  if (num_remaining_args > 0) {
    if (OPS_GetDoubleInput(&num_remaining_args, b_data) < 0) {
      opserr << "WARNING: invalid double b_data\n";
      return 0;
    }
  }

  theElement = new ShellDKGT(iData[0], iData[1], iData[2], iData[3],
                             *theSection, b_data[0], b_data[1], b_data[2]);

  return theElement;
}
