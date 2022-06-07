

#include <g3_api.h>


#include <SRC/element/frictionBearing/frictionModel/Coulomb.h>
void *OPS_Coulomb(void)
{
  // pointer to a friction model that will be returned
  FrictionModel *theFrnMdl = 0;

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: frictionModel Coulomb tag mu\n";
    return 0;
  }

  int tag[1];
  double dData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid tag for frictionModel Coulomb\n";
    return 0;
  }
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for frictionModel Coulomb " << tag[0] << endln;
    return 0;
  }

  // parsing was successful, allocate the friction model
  theFrnMdl = new Coulomb(tag[0], dData[0]);
  if (theFrnMdl == 0) {
    opserr << "WARNING could not create frictionModel of type Coulomb\n";
    return 0;
  }

  return theFrnMdl;
}
