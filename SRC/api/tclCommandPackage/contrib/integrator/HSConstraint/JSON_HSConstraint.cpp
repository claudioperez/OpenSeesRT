

#include <g3_api.h>


#include <SRC/analysis/integrator/HSConstraint.h>
void *OPS_HSConstraint() {
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 1) {
    opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> "
              "<u_ref> \n";
    return 0;
  }
  if (numdata > 4)
    numdata = 4;

  double data[4];
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING integrator HSConstraint invalid double inputs\n";
    return 0;
  }
  double arcLength = data[0];
  double psi_u = data[1];
  double psi_f = data[2];
  double u_ref = data[3];

  switch (numdata) {
  case 1:
    return new HSConstraint(arcLength);
  case 2:
    return new HSConstraint(arcLength, psi_u);
  case 3:
    return new HSConstraint(arcLength, psi_u, psi_f);
  case 4:
    return new HSConstraint(arcLength, psi_u, psi_f, u_ref);
  }

  return 0;
}
