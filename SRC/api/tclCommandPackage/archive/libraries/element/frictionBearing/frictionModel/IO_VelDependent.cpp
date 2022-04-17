

#include <g3_api.h>


#include <SRC/element/frictionBearing/frictionModel/VelDependent.h>
void *OPS_VelDependent(void)
{
  // pointer to a friction model that will be returned
  FrictionModel *theFrnMdl = 0;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: frictionModel VelDependent tag muSlow muFast transRate\n";
    return 0;
  }

  int tag[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid tag for frictionModel VelDependent\n";
    return 0;
  }
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for frictionModel VelDependent " << tag[0] << endln;
    return 0;
  }

  // parsing was successful, allocate the friction model
  theFrnMdl = new VelDependent(tag[0], dData[0], dData[1], dData[2]);
  if (theFrnMdl == 0) {
    opserr << "WARNING could not create frictionModel of type VelDependent\n";
    return 0;
  }

  return theFrnMdl;
}
