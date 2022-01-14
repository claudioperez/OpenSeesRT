

#include <g3_api.h>


#include <SRC/element/frictionBearing/frictionModel/VelDepMultiLinear.h>
void *OPS_VelDepMultiLinear(void)
{
  // pointer to a friction model that will be returned
  FrictionModel *theFrnMdl = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 7) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: frictionModel VelDepMultiLinear tag ";
    opserr << "-vel velocityPoints -frn frictionPoints  ";
    opserr << "(with at least two friction-velocity points)";
    return 0;
  }

  int tag[1];
  double velData[64];
  double frnData[64];
  const char *paraStr;
  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid tag for frictionModel VelDepMultiLinear\n";
    return 0;
  }

  // get velocity data points
  numData = (argc - 3) / 2;
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-vel") == 0) {
    if (OPS_GetDoubleInput(&numData, velData) != 0) {
      opserr << "WARNING invalid velocityPoints\n";
      opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -vel but got " << paraStr << endln;
    opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
    return 0;
  }
  Vector velPts(velData, numData);

  // get friction data points
  paraStr = OPS_GetString();
  if (strcmp(paraStr, "-frn") == 0) {
    if (OPS_GetDoubleInput(&numData, frnData) != 0) {
      opserr << "WARNING invalid frictionPoints\n";
      opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
      return 0;
    }
  } else {
    opserr << "WARNING expecting -frn but got " << paraStr << endln;
    opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
    return 0;
  }
  Vector frnPts(frnData, numData);

  // parsing was successful, allocate the friction model
  theFrnMdl = new VelDepMultiLinear(tag[0], velPts, frnPts);
  if (theFrnMdl == 0) {
    opserr
        << "WARNING could not create frictionModel of type VelDepMultiLinear\n";
    return 0;
  }

  return theFrnMdl;
}
