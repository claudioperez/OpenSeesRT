

#include <g3_api.h>


#include <SRC/analysis/integrator/BackwardEuler.h>
void *OPS_BackwardEuler() {
  int optn = 0;
  if (OPS_GetNumRemainingInputArgs() > 0) {
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &optn) < 0) {
      opserr << "WARNING integrator BackwardEuler <option> - undefined option "
                "specified\n";
      return 0;
    }
  }
  return new BackwardEuler(optn);
}
