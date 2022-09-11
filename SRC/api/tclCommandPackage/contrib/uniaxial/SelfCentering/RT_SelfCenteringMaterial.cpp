#include "SelfCenteringMaterial.h"
#include <InputAPI.h>

UniaxialMaterial*
G3Parse_newSelfCenteringMaterial(ClientData* cd, Tcl_Interp* interp, int argc, G3_Char** argv)
{
  UniaxialMaterial* theMaterial;
     //if (strcmp(argv[1], "SelfCentering") == 0) {
      if (argc < 7) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr
            << "Want: uniaxialMaterial SelfCentering tag? k1? k2? ActF? beta? "
               "<SlipDef? BearDef? rBear?>"
            << endln;
        return nullptr;
      }

      int tag;
      double k1, k2, ActF, beta, rBear, SlipDef, BearDef;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial SelfCentering tag" << endln;
        return nullptr;
      }

      if (Tcl_GetDouble(interp, argv[3], &k1) != TCL_OK) {
        opserr << "WARNING invalid k1\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return nullptr;
      }

      if (Tcl_GetDouble(interp, argv[4], &k2) != TCL_OK) {
        opserr << "WARNING invalid k2\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return nullptr;
      }

      if (Tcl_GetDouble(interp, argv[5], &ActF) != TCL_OK) {
        opserr << "WARNING invalid ActF\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return nullptr;
      }

      if (Tcl_GetDouble(interp, argv[6], &beta) != TCL_OK) {
        opserr << "WARNING invalid beta\n";
        opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
        return nullptr;
      }

      if (argc == 8) {
        if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
          opserr << "WARNING invalid SlipDef\n";
          opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
          return nullptr;
        }
        // Parsing was successful, allocate the material
        theMaterial =
            new SelfCenteringMaterial(tag, k1, k2, ActF, beta, SlipDef, 0, 0);
      }

      else if (argc > 8) {
        if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
          opserr << "WARNING invalid SlipDef\n";
          opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[8], &BearDef) != TCL_OK) {
          opserr << "WARNING invalid BearDef\n";
          opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[9], &rBear) != TCL_OK) {
          opserr << "WARNING invalid rBear\n";
          opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
          return nullptr;
        }
        // Parsing was successful, allocate the material
        theMaterial = new SelfCenteringMaterial(tag, k1, k2, ActF, beta,
                                                SlipDef, BearDef, rBear);
      }

      else {
        // Parsing was successful, allocate the material
        theMaterial =
            new SelfCenteringMaterial(tag, k1, k2, ActF, beta, 0, 0, 0);
      }
  return theMaterial;
    ///}
}
