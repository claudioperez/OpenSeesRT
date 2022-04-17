

#include <g3_api.h>


#include <SRC/element/frictionBearing/frictionModel/VelPressureDep.h>
void *OPS_VelPressureDep(void)
{
  // pointer to a friction model that will be returned
  FrictionModel *theFrnMdl = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: frictionModel VelPressureDep tag muSlow muFast0 A deltaMu "
              "alpha transRate\n";
    return 0;
  }

  int tag[1];
  double dData[6];
  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid tag for frictionModel VelPressureDep\n";
    return 0;
  }
  numData = 6;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for frictionModel VelPressureDep " << tag[0]
           << endln;
    return 0;
  }

  // parsing was successful, allocate the friction model
  theFrnMdl = new VelPressureDep(tag[0], dData[0], dData[1], dData[2], dData[3],
                                 dData[4], dData[5]);
  if (theFrnMdl == 0) {
    opserr << "WARNING could not create frictionModel of type VelPressureDep\n";
    return 0;
  }

  return theFrnMdl;
}
