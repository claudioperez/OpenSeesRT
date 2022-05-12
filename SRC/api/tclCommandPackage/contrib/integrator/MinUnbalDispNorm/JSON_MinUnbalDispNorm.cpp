

#include <g3_api.h>


#include <SRC/analysis/integrator/MinUnbalDispNorm.h>
void *OPS_MinUnbalDispNorm() {
  double lambda11, minlambda, maxlambda;
  int numIter;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j "
              "maxLambda1j>\n";
    return 0;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &lambda11) < 0) {
    opserr << "WARNING integrator MinUnbalDispNorm invalid lambda11\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() >= 3) {
    if (OPS_GetIntInput(&numdata, &numIter) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid numIter\n";
      return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &minlambda) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid minlambda\n";
      return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &maxlambda) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid maxlambda\n";
      return 0;
    }
  } else {
    minlambda = lambda11;
    maxlambda = lambda11;
    numIter = 1;
  }

  int signFirstStepMethod = SIGN_LAST_STEP;
  if (OPS_GetNumRemainingInputArgs() > 0) {
    const char *flag = OPS_GetString();
    if ((strcmp(flag, "-determinant") == 0) || (strcmp(flag, "-det") == 0)) {
      signFirstStepMethod = CHANGE_DETERMINANT;
    }
  }

  return new MinUnbalDispNorm(lambda11, numIter, minlambda, maxlambda,
                              signFirstStepMethod);
}
