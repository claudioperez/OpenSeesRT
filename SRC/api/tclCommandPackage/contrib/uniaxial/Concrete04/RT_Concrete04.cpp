#include <G3Parse.h>
#include <Concrete04.h>

UniaxialMaterial* G3Parse_newUniaxialConcrete04(G3_Runtime* rt, int argc, G3_Char** argv)
/* if (strcmp(argv[1], "Concrete04") == 0) */
{
    UniaxialMaterial *theMaterial = 0;
    int tag;
    double fpc, epsc0, ft, epscu, Ec0, etu, beta;

    if (argc != 10 && argc != 9 && argc != 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? "
                "Ec0? <ft? etu? <beta?> >"
             << endln;
      return nullptr;
    }


    if (G3Parse_getInt(rt, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Concrete04 tag" << endln;
      return nullptr;
    }

    // Read required Concrete04 material parameters

    if (G3Parse_getDouble(rt, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      opserr << "Concrete04 material: " << tag << endln;
      return nullptr;
    }

    if (G3Parse_getDouble(rt, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return nullptr;
    }

    if (G3Parse_getDouble(rt, argv[5], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "Concrete04 material: " << tag << endln;
      return nullptr;
    }

    if (G3Parse_getDouble(rt, argv[6], &Ec0) != TCL_OK) {
      opserr << "WARNING invalid Ec0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return nullptr;
    }
    if (argc == 9 || argc == 10) {
      if (G3Parse_getDouble(rt, argv[7], &ft) != TCL_OK) {
        opserr << "WARNING invalid ft\n";
        opserr << "Concrete04 material: " << tag << endln;
        return nullptr;
      }
      if (G3Parse_getDouble(rt, argv[8], &etu) != TCL_OK) {
        opserr << "WARNING invalid etu\n";
        opserr << "Concrete04 material: " << tag << endln;
        return nullptr;
      }
    }
    if (argc == 10) {
      if (G3Parse_getDouble(rt, argv[9], &beta) != TCL_OK) {
        opserr << "WARNING invalid beta\n";
        opserr << "Concrete04 material: " << tag << endln;
        return nullptr;
      }
    }

    // Parsing was successful, allocate the material
    if (argc == 10) {
      theMaterial =
          new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu, beta);
    } else if (argc == 9) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu);
    } else if (argc == 7) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0);
    }
    //return G3_addUniaxialMaterial(G3_getRuntime(interp), theMaterial);
    return theMaterial;
}

