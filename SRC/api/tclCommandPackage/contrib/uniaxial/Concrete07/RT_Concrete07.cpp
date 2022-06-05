
#include <G3Parse.h>
#include <Concrete07.h>

UniaxialMaterial*
G3Parse_newUniaxialConcrete07(G3_Runtime* rt, int argc, G3_Char** argv)
/* else if (strcmp(argv[1], "Concrete07") == 0) */ 
{
    UniaxialMaterial *theMaterial = 0;
      // Check to see if there are enough arquements
      if (argc < 11) {
        opserr << "WARNING: Insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: uniaxialMaterial Concrete07 tag? fpc? epsc0? Ec? fpt? "
                  "epst0? xcrp? xcrn? r?\n";
        return nullptr;
      }

      int tag;

      if (G3Parse_getInt(rt, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING: Invalid uniaxial Concrete07 tag\n";
        return nullptr;
      }

      // Read in the faluves required for the model
      double fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r;

      if (G3Parse_getDouble(rt, argv[3], &fpc) != TCL_OK) {
        opserr << "WARNING: Invalid peak compression stress\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[4], &epsc0) != TCL_OK) {
        opserr << "WARNING: Invalid peak compression strain\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[5], &Ec) != TCL_OK) {
        opserr << "WARNING: Invalid Young's Modulus\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[6], &fpt) != TCL_OK) {
        opserr << "WARNING: Invalid peak tension stress\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[7], &epst0) != TCL_OK) {
        opserr << "WARNING: Invalid peak tension strain\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[8], &xcrp) != TCL_OK) {
        opserr
            << "WARNING: Invalid critical nondimensional strain in tension\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[9], &xcrn) != TCL_OK) {
        opserr << "WARNING: Invalid critical nondimensional strain in "
                  "compression\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[10], &r) != TCL_OK) {
        opserr << "WARNING: Invalid value for r\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      }

      // Parsing was successful, allocate the material
      theMaterial =
          new Concrete07(tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);

    //return G3_addUniaxialMaterial(G3_getRuntime(interp), theMaterial);
    return theMaterial;
}   
