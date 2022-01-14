

#include <g3_api.h>


#include <SRC/element/componentElement/ComponentElement2d.h>
void *OPS_ComponentElement2d(void)
{
  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr << "Invalid #args,  want: element CompositeElement tag iNode jNode "
              "A E I crdTag hinge1 hinge2 \n";
    return 0;
  }

  int iData[6];
  double dData[3];
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids ints" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids double" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids second set ints" << endln;
    return 0;
  }

  double mass = 0.0;
  int cMass = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if (type == "-rho") {
      int numData = 1;
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0)
          return 0;
      }
    } else if (type == "-cMass") {
      cMass = 1;
    }
  }

  CrdTransf *theTrans = OPS_getCrdTransf(iData[3]);

  UniaxialMaterial *end1 = OPS_getUniaxialMaterial(iData[4]);
  UniaxialMaterial *end2 = OPS_getUniaxialMaterial(iData[5]);

  // Parsing was successful, allocate the material
  theElement =
      new ComponentElement2d(iData[0], dData[0], dData[1], dData[2], iData[1],
                             iData[2], *theTrans, end1, end2, mass, cMass);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type ComponentElement2d\n";
    return 0;
  }

  return theElement;
}
