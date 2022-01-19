

#include <g3_api.h>


#include <SRC/element/truss/TrussSection.h>
OPS_Export void *OPS_TrussSectionElement()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element TrussSection $tag $iNode $jNode "
              "$sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }

  int iData[4];
  double rho = 0.0;
  int ndm = OPS_GetNDM();
  int doRayleigh = 0; // by default rayleigh not done
  int cMass = 0;      // by default use lumped mass matrix

  int numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, sectTag) in element "
              "TrussSection "
           << endln;
    return 0;
  }

  SectionForceDeformation *theSection =
      OPS_getSectionForceDeformation(iData[3]);

  if (theSection == 0) {
    opserr << "WARNING: Invalid section not found element TrussSection "
           << iData[0] << " $iNode $jNode " << iData[3]
           << " <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }

  numRemainingArgs -= 4;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();

    if (strcmp(argvS, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
        opserr << "WARNING Invalid rho in element TrussSection " << iData[0]
               << " $iNode $jNode $secTag <-rho $rho> <-cMass $flag> "
                  "<-doRayleigh $flag>\n";
        return 0;
      }
    } else if (strcmp(argvS, "-cMass") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &cMass) != 0) {
        opserr << "WARNING: Invalid cMass in element TrussSection " << iData[0]
               << " $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> "
                  "<-doRayleigh $flag>\n";
        return 0;
      }
    } else if (strcmp(argvS, "-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
        opserr << "WARNING: Invalid doRayleigh in element TrussSection "
               << iData[0]
               << " $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> "
                  "<-doRayleigh $flag>\n";
        return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS
             << "  in: element TrussSection " << iData[0]
             << " $iNode $jNode $secTag <-rho $rho> <-cMass $flag> "
                "<-doRayleigh $flag>\n";
      return 0;
    }
    numRemainingArgs -= 2;
  }

  // now create the TrussSection
  theElement = new TrussSection(iData[0], ndm, iData[1], iData[2], *theSection,
                                rho, doRayleigh, cMass);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element TrussSection " << iData[0]
           << " $iNode $jNode $secTag <-rho $rho> <-cMass $flag> <-doRayleigh "
              "$flag>\n";
  }

  return theElement;
}
