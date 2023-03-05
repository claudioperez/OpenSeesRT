

#include <g3_api.h>


#include <SRC/analysis/integrator/DisplacementControl.h>
void *OPS_DisplacementControlIntegrator() {
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "insufficient arguments for DisplacementControl\n";
    return 0;
  }

  // node, dof
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING failed to read node tag and ndf\n";
    return 0;
  }

  double incr;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &incr) < 0) {
    opserr << "WARNING failed to read incr\n";
    return 0;
  }

  // numIter,dumin,dumax
  int numIter = 1;
  int formTangent = 0;
  double data[2] = {incr, incr};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING failed to read numIter\n";
      return 0;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
      opserr << "WARNING failed to read dumin and dumax\n";
      return 0;
    }
  }

  if (OPS_GetNumRemainingInputArgs() == 1) {
    std::string type = OPS_GetString();
    if (type == "-initial" || type == "-Initial") {
      formTangent = 1;
    }
  }

  // check node
  Domain *theDomain = OPS_GetDomain();
  Node *theNode = theDomain->getNode(iData[0]);
  if (theNode == 0) {
    opserr << "WARNING integrator DisplacementControl node dof dU : Node does "
              "not exist\n";
    return 0;
  }

  int numDOF = theNode->getNumberDOF();
  if (iData[1] <= 0 || iData[1] > numDOF) {
    opserr << "WARNING integrator DisplacementControl node dof dU : invalid "
              "dof given\n";
    return 0;
  }

  return new DisplacementControl(iData[0], iData[1] - 1, incr, theDomain,
                                 numIter, data[0], data[1], formTangent);
}
