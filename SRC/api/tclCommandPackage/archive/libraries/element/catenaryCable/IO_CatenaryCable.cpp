

#include <g3_api.h>


#include <SRC/element/catenaryCable/CatenaryCable.h>
OPS_Export void *OPS_CatenaryCableElement()
{

  if (num_CatenaryCableElement == 0) {
    num_CatenaryCableElement++;
    opserr << "CatenaryCableElement element - Written: P. Ibanez and J. A. "
              "Abell (UANDES). www.joseabell.com.\n";
  }

  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  //(int tag, int node1, int node2, double weight, double E, double A, double
  //L0, double alpha, double temperature_change = 100., double rho)
  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element CatenaryCable $tag $iNode $jNode "
              "$weight $E $A $L0 $alpha $temperature_change $rho $errorTol "
              "$Nsubsteps\n";
    return 0;
  }

  if (numRemainingArgs != 13) {
    opserr << "Got " << numRemainingArgs << " args. Expected 13\n";
    opserr << "Invalid Args want: element CatenaryCable $tag $iNode $jNode "
              "$weight $E $A $L0 $alpha $temperature_change $rho $errorTol "
              "$Nsubsteps $massType\n";
    return 0; // it's a CatenaryCableSection
  }

  int iData[3];
  double dData[8];

  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer (tag, iNode, "
              "jNode) in element CatenaryCable "
           << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING:  element CatenaryCable - invalid double data. Check "
              "$weight $E $A $L0 $alpha $temperature_change $rho $errorTol "
              "$Nsubsteps $massType\n";
    return 0;
  }

  numData = 1;
  int Nsubsteps = 0;
  if (OPS_GetInt(&numData, &Nsubsteps) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer $Nsubsteps in "
              "element CatenaryCable "
           << endln;
    return 0;
  }

  int massType = 0;
  if (OPS_GetInt(&numData, &massType) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer $massType in "
              "element CatenaryCable "
           << endln;
    return 0;
  }

  // now create the CatenaryCable
  theElement = new CatenaryCable(
      iData[0], iData[1], iData[2], dData[0], dData[1], dData[2], dData[3],
      dData[4], dData[5], dData[6], dData[7], Nsubsteps, massType);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CatenaryCable " << iData[0]
           << " $iNode $jNode ...\n";
  }

  return theElement;
}
