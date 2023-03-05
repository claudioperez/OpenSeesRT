

#include <g3_api.h>


#include <SRC/analysis/integrator/HHTGeneralizedExplicit_TP.h>
void *OPS_HHTGeneralizedExplicit_TP(void) {
  // pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want "
              "HHTGeneralizedExplicit_TP $rhoB $alphaF\n";
    opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta "
              "$gamma\n";
    return 0;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralizedExplicit_TP $rhoB "
              "$alphaF\n";
    opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta "
              "$gamma\n";
    return 0;
  }

  if (argc == 2)
    theIntegrator = new HHTGeneralizedExplicit_TP(dData[0], dData[1]);
  else if (argc == 4)
    theIntegrator =
        new HHTGeneralizedExplicit_TP(dData[0], dData[1], dData[2], dData[3]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating HHTGeneralizedExplicit_TP "
              "integrator\n";

  return theIntegrator;
}
