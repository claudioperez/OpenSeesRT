

#include <g3_api.h>


#include <SRC/element/truss/CorotTruss2.h>
void *OPS_CorotTruss2(void)
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 7) {
    opserr << "Invalid Args want: element CorotTruss2 $tag $iNode $jNode "
              "$auxN1 $auxN2 $A $matTag <-rho $rho>\n";
    return 0;
  }

  int iData[5];
  double A = 0.0;
  double rho = 0.0;
  int matTag = 0;
  int ndm = OPS_GetNDM();

  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, auxN1, auxN2) in "
              "element CorotTruss2 "
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element CorotTruss2 " << iData[0]
           << " $iNode $jNode $A $matTag <-rho $rho> <-rayleig $flagh>\n";
    return 0;
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag) != 0) {
    opserr << "WARNING: Invalid matTag: element CorotTruss2 " << iData[0]
           << " $iNode $jNode $A $matTag <-rho $rho> <-rayleig $flagh>\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);

  if (theUniaxialMaterial == 0) {
    opserr << "WARNING: Invalid material not found element CorotTruss2 "
           << iData[0] << " $iNode $jNode $auxN1 $auxN2 $A " << matTag
           << " <-rho $rho> <-rayleigh $flagh>\n";
    return 0;
  }

  numRemainingArgs -= 7;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();

    if (strcmp(argvS, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
        opserr << "WARNING Invalid rho in element CorotTruss2 " << iData[0]
               << " $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> "
                  "<-rayleigh $flagh>\n";
        return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS
             << "  in: element CorotTruss2 " << iData[0]
             << " $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> "
                "<-rayleigh $flagh>\n";
      return 0;
    }
    numRemainingArgs -= 2;
  }

  // now create the ReinforcedConcretePlaneStress
  theElement = new CorotTruss2(iData[0], ndm, iData[1], iData[2], iData[3],
                               iData[4], *theUniaxialMaterial, A, rho);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CorotTruss2 " << iData[0]
           << " $iNode $jNode $A $matTag <-rho $rho> \n";
  }

  return theElement;
}
