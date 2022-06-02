#include <G3Parse.h>
#include <BoucWenMaterial.h>

UniaxialMaterial*
G3Parse_newUniaxialBoucWen(G3_Builder* rt, int argc, G3_Char** argv)
{
  // else if (strcmp(argv[1], "BoucWen") == 0) {
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?"
             << endln << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
      return nullptr;
    }

    int tag;
    double alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta;

    if (G3Parse_getInt(rt, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial BoucWen tag" << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[3], &alpha) != TCL_OK) {
      opserr << "WARNING invalid alpha\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[4], &ko) != TCL_OK) {
      opserr << "WARNING invalid ko\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[5], &n) != TCL_OK) {
      opserr << "WARNING invalid n\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[6], &gamma) != TCL_OK) {
      opserr << "WARNING invalid gamma\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[7], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[8], &Ao) != TCL_OK) {
      opserr << "WARNING invalid Ao\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[9], &deltaA) != TCL_OK) {
      opserr << "WARNING invalid deltaA\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[10], &deltaNu) != TCL_OK) {
      opserr << "WARNING invalid deltaNu\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }
    if (G3Parse_getDouble(rt, argv[11], &deltaEta) != TCL_OK) {
      opserr << "WARNING invalid deltaEta\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return nullptr;
    }

    // Check if the user has given a tolerance for the Newton scheme
    double tolerance = 1.0e-8;
    if (argc > 12) {
      if (G3Parse_getDouble(rt, argv[12], &tolerance) != TCL_OK) {
        opserr << "WARNING invalid tolerance\n";
        opserr << "uniaxialMaterial BoucWen: " << tolerance << endln;
        return nullptr;
      }
    }

    // Check if the user has given a maxNumIter for the Newton scheme
    int maxNumIter = 20;
    if (argc > 13) {
      if (G3Parse_getInt(rt, argv[13], &maxNumIter) != TCL_OK) {
        opserr << "WARNING invalid maxNumIter\n";
        opserr << "uniaxialMaterial BoucWen: " << maxNumIter << endln;
        return nullptr;
      }
    }

    // Parsing was successful, allocate the material
    return
        new BoucWenMaterial(tag, alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu,
                            deltaEta, tolerance, maxNumIter);
}
