

#include <g3_api.h>


#include <SRC/element/adapter/ActuatorCorot.h>
void *OPS_ActuatorCorot()
{
  // check the number of arguments is correct
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element actuator eleTag iNode jNode EA ipPort <-ssl> "
              "<-udp> <-doRayleigh> <-rho rho>\n";
    return 0;
  }

  int ndm = OPS_GetNDM();

  // get the id and end nodes
  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid actuator int inputs" << endln;
    return 0;
  }
  int tag = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];

  // stiffness
  double EA;
  numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &EA) < 0) {
    opserr << "WARNING invalid actuator EA" << endln;
    return 0;
  }

  // ipPort
  int ipPort;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &ipPort) < 0) {
    opserr << "WARNING invalid actuator ipPort" << endln;
    return 0;
  }

  // options
  int ssl = 0, udp = 0;
  int doRayleigh = 0;
  double rho = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *flag = OPS_GetString();
    if (strcmp(flag, "-ssl") == 0) {
      ssl = 1;
      udp = 0;
    } else if (strcmp(flag, "-udp") == 0) {
      udp = 1;
      ssl = 0;
    } else if (strcmp(flag, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(flag, "-rho") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        numdata = 1;
        if (OPS_GetDoubleInput(&numdata, &rho) < 0) {
          opserr << "WARNING invalid rho\n";
          opserr << "actuator element: " << tag << endln;
          return 0;
        }
      }
    }
  }

  // now create the actuator and add it to the Domain
  return new ActuatorCorot(tag, ndm, iNode, jNode, EA, ipPort, ssl, udp,
                           doRayleigh, rho);
}
