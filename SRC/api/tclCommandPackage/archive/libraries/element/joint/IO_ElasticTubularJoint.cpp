

#include <g3_api.h>


#include <SRC/element/joint/ElasticTubularJoint.h>
void *OPS_ElasticTubularJoint(void)
{

  if (numElasticTubularJoint == 0) {
    numElasticTubularJoint++;
    opserr << "ElasticTubularJoint element - Written by Kia & Alanjari\n";
  }

  // get the id and end nodes
  int iData[3];
  double dData[6];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
    opserr << "\n WARNING invalid ElasticTubularJoint Tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "\n WARNING invalid iNode for ElasticTubularJoint " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
    opserr << "\n WARNING invalid jNode for ElasticTubularJoint " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
    opserr << "\n WARNING invalid  brace diameter for ElasticTubularJoint "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "\n WARNING invalid  brace_angle for ElasticTubularJoint "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
    opserr << "\n WARNING invalid E  for ElasticTubularJoint " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
    opserr << "\n WARNING invalid  chord diameter for ElasticTubularJoint "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[4]) != 0) {
    opserr << "\n WARNING invalid  chord thickness for ElasticTubularJoint "
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[5]) != 0) {
    opserr << "\n WARNING invalid  chord angle for ElasticTubularJoint "
           << iData[0] << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  // Element *theTubularJoint = new ElasticTubularJoint(iData[0], iData[1],
  // iData[2],dData[3] ,dData[2] , dData[5] , dData[0] ,dData[1] ,dData[4]);

  Element *theTubularJoint =
      new ElasticTubularJoint(iData[0], iData[1], iData[2], dData[0], dData[1],
                              dData[2], dData[3], dData[4], dData[5]);

  if (theTubularJoint == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << iData[0]
           << endln;
    return 0;
  }

  return theTubularJoint;
}
