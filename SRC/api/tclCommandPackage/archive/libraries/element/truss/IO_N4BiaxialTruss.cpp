

#include <g3_api.h>


#include <SRC/element/truss/N4BiaxialTruss.h>
void *OPS_N4BiaxialTruss(void)
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 7) {
    opserr << "Invalid Args want: element N4BiaxialTruss $tag $i1Node $j1Node "
              "$iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
    return 0;
  }

  int iData[5]; // tag, iNode, jNode, iGNode, jGNode
  double A = 0.0;
  double rho = 0.0;
  int matTag1 = 0;
  int matTag2 = 0;
  int doRayleigh = 0;
  int ndm = OPS_GetNDM();

  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, iGNode, jGNode) in "
              "element N4BiaxialTruss "
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element N4BiaxialTruss " << iData[0]
           << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
              "<-doRayleigh $flag>\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag1) != 0) {
    opserr << "WARNING: Invalid matTag1: element N4BiaxialTruss " << iData[0]
           << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
              "<-doRayleigh $flag>\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial_1 = OPS_GetUniaxialMaterial(matTag1);
  if (theUniaxialMaterial_1 == 0) {
    opserr << "WARNING: Invalid material not found element N4BiaxialTruss "
           << iData[0] << " $mattag1: " << matTag1 << " \n";
    return 0;
  }

  numRemainingArgs -= 6;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();

    if (strcmp(argvS, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
        opserr << "WARNING Invalid rho in element N4BiaxialTruss " << iData[0]
               << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
                  "<-doRayleigh $flag>\n";
        return 0;
      }
    } else if (strcmp(argvS, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
        opserr << "WARNING: Invalid doRayleigh in element N4BiaxialTruss "
               << iData[0]
               << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
                  "<-doRayleigh $flag>\n";
        return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS
             << "  in: element N4BiaxialTruss " << iData[0]
             << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
                "<-doRayleigh $flag>\n";
      return 0;
    }
    numRemainingArgs -= 2;
  }

  // now create the ReinforcedConcretePlaneStress
  theElement =
      new N4BiaxialTruss(iData[0], ndm, iData[1], iData[2], iData[3], iData[4],
                         *theUniaxialMaterial_1, A, rho, doRayleigh);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element N4BiaxialTruss " << iData[0]
           << " $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> "
              "<-doRayleigh $flag>\n";
  }

  return theElement;
}
