

#include <g3_api.h>


#include <SRC/element/frictionBearing/frictionModel/VelNormalFrcDep.h>
void *OPS_VelNormalFrcDep(void)
{
  // pointer to a friction model that will be returned
  FrictionModel *theFrnMdl = 0;

  if (OPS_GetNumRemainingInputArgs() < 9) {
    opserr << "WARNING invalid number of arguments\n";
    opserr << "Want: frictionModel VelNormalFrcDep tag aSlow nSlow aFast nFast "
              "alpha0 alpha1 alpha2 maxMuFact\n";
    return 0;
  }

  int tag[1];
  double dData[8];
  int numData = 1;
  if (OPS_GetIntInput(&numData, tag) != 0) {
    opserr << "WARNING invalid tag for frictionModel VelNormalFrcDep\n";
    return 0;
  }
  numData = 8;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for frictionModel VelNormalFrcDep " << tag[0]
           << endln;
    return 0;
  }

  // parsing was successful, allocate the friction model
  theFrnMdl =
      new VelNormalFrcDep(tag[0], dData[0], dData[1], dData[2], dData[3],
                          dData[4], dData[5], dData[6], dData[7]);
  if (theFrnMdl == 0) {
    opserr
        << "WARNING could not create frictionModel of type VelNormalFrcDep\n";
    return 0;
  }

  return theFrnMdl;
}
