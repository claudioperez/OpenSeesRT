

#include <g3_api.h>


#include <SRC/analysis/integrator/TRBDF2.h>
void *OPS_TRBDF2() { return new TRBDF2(); }

TRBDF2::TRBDF2()
    : TransientIntegrator(INTEGRATOR_TAGS_TRBDF2), step(0), dt(0.0), c1(0.0),
      c2(0.0), c3(0.0), Utm1(0), Utm1dot(0), Ut(0), Utdot(0), Utdotdot(0), U(0),
      Udot(0), Udotdot(0) {}
TRBDF2::~TRBDF2() {
  // clean up the memory created
  if (Utm1 != 0)
    delete Utm1;
  if (Utm1dot != 0)
    delete Utm1dot;

  if (Ut != 0)
    delete Ut;
  if (Utdot != 0)
    delete Utdot;
  if (Utdotdot != 0)
    delete Utdotdot;

  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
}
