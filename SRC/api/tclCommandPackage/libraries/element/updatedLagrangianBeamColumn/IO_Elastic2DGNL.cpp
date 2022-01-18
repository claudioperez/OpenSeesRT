

#include <g3_api.h>


#include <SRC/element/updatedLagrangianBeamColumn/Elastic2DGNL.h>
void *OPS_Elastic2DGNL()
{
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dGNL int tag, int Nd1, int Nd2, double A, "
              "double E, double Iz, <int linear>\n";

    return 0;
  }

  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid Elastic2dGNL int inputs" << endln;
    return 0;
  }

  int tag = idata[0];
  int ndI = idata[1];
  int ndJ = idata[2];

  double data[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid Elastic2dGNL double inputs" << endln;
    return 0;
  }

  double A = data[0];
  double E = data[1];
  double I = data[2];

  bool linear = false;

  if (OPS_GetNumRemainingInputArgs() > 0) {
    numdata = 1;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
      opserr << "WARNING invalid Elastic2dGNL int inputs" << endln;
      return 0;
    }
    if (idata[0] == 1)
      linear = true;
  }

  return new Elastic2dGNL(tag, A, E, I, ndI, ndJ, linear); //, false, massDens);
}
