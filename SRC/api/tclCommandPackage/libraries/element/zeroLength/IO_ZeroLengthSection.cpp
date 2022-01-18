

#include <g3_api.h>


#include <SRC/element/zeroLength/ZeroLengthSection.h>
void *OPS_ZeroLengthSection()
{
  int ndm = OPS_GetNDM();

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "insufficient arguments for ZeroLengthSection\n";
    return 0;
  }

  // get eleTag,iNode,jNode,secTag
  int iData[4];
  int numData = 4;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING: invalid integer inputs\n";
    return 0;
  }

  // options
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;
  double *x_ptr = &x(0), *y_ptr = &y(0);
  int doRayleighDamping = 1;
  while (OPS_GetNumRemainingInputArgs() > 1) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-orient") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 5) {
        numData = 3;
        if (OPS_GetDoubleInput(&numData, x_ptr) < 0) {
          opserr << "WARNING: invalid double inputs\n";
          return 0;
        }
        if (OPS_GetDoubleInput(&numData, y_ptr) < 0) {
          opserr << "WARNING: invalid double inputs\n";
          return 0;
        }
      }
    } else if (strcmp(type, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &doRayleighDamping) < 0) {
        opserr << "WARNING: invalid integer inputs\n";
        return 0;
      }
    }
  }

  // get section
  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[3]);
  if (theSection == 0) {
    opserr << "zeroLengthSection -- no section with tag " << iData[0]
           << " exists in Domain\n";
    return 0;
  }

  return new ZeroLengthSection(iData[0], ndm, iData[1], iData[2], x, y,
                               *theSection, doRayleighDamping);
}
