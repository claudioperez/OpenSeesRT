

#include <g3_api.h>


#include <SRC/analysis/integrator/Newmark1.h>
void *OPS_Newmark1() {
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata != 2 && numdata != 6) {
    opserr << "WARNING integrator Newmark1 gamma beta <alphaM> <betaKcurrent> "
              "<betaKi> <betaKlastCommitted>\n";
    return 0;
  }

  double data[6] = {0, 0, 0, 0, 0, 0};
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING integrator Newmark1 invalid double inputs\n";
    return 0;
  }

  double gamma = data[0];
  double beta = data[1];
  double alphaM = data[2], betaK = data[3], betaKi = data[4], betaKc = data[5];

  if (numdata == 2)
    return new Newmark1(gamma, beta);
  else
    return new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc);
}
