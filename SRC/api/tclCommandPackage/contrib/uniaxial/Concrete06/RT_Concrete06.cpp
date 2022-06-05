
#include <G3Parse.h>
#include <Concrete06.h>

/*
static void printCommand(int argc, G3_Char **argv) {
  opserr << "Input command: ";
  for (int i = 0; i < argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
}
*/

UniaxialMaterial* G3Parse_newUniaxialConcrete06(G3_Runtime* rt, int argc, G3_Char** argv)
{
  UniaxialMaterial *theMaterial = 0;
      if (argc < 12) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr
            << "Want: uniaxialMaterial Concrete06 tag? fc? eo? r? k? alphaC? "
               "fcr? ecr? b? alphaT?"
            << endln;
        return nullptr;
      }

      int tag;

      if (G3Parse_getInt(rt, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Concrete06 tag" << endln;
        return nullptr;
      }

      // Read required Concrete01 material parameters
      double fc, eo, r, k, fcr, ecr, b, alphaC, alphaT;

      if (G3Parse_getDouble(rt, argv[3], &fc) != TCL_OK) {
        opserr << "WARNING invalid fc\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[4], &eo) != TCL_OK) {
        opserr << "WARNING invalid eo\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[5], &r) != TCL_OK) {
        opserr << "WARNING invalid r\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[6], &k) != TCL_OK) {
        opserr << "WARNING invalid k\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[7], &alphaC) != TCL_OK) {
        opserr << "WARNING invalid alphaC\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[8], &fcr) != TCL_OK) {
        opserr << "WARNING invalid fcr\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[9], &ecr) != TCL_OK) {
        opserr << "WARNING invalid ecr\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[10], &b) != TCL_OK) {
        opserr << "WARNING invalid b\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      if (G3Parse_getDouble(rt, argv[11], &alphaT) != TCL_OK) {
        opserr << "WARNING invalid alphaT\n";
        opserr << "Concrete06 material: " << tag << endln;
        return nullptr;
      }

      // Parsing was successful, allocate the material
      theMaterial =
          new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);

    //return G3_addUniaxialMaterial(G3_getRuntime(interp), theMaterial);
    return theMaterial;
}


