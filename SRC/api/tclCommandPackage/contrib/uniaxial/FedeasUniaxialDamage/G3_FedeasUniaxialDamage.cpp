#include <stdio.h>
#include <G3Parse.h>
#include <DegradingUniaxialWrapper.hh>

#define WRAPPER_CMD "FedeasUniaxialDamage"
// #define WRAPPER_CMD "FedeasDamage"

UniaxialMaterial*
G3Parse_newFedeasUniaxialDamage(G3_Runtime* rt, int argc, TCL_Char **argv)
{
  // Pointer to a uniaxial material that will be returned
  DegradingUniaxialWrapper *theMaterial = 0;
  UniaxialMaterial *theWrappedMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain =  1.0e16;
  int tags[2];

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-damage $damageTag>"
           << endln;
    return nullptr;
  }

  // Get wrapper tag
  if (G3Parse_getInt(rt, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }
  // Get base tag
  if (G3Parse_getInt(rt, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }

  // Get base material
  theWrappedMaterial = G3_getUniaxialMaterialInstance(rt, tags[1]);
  if (theWrappedMaterial == 0) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag" WRAPPER_CMD " tag: "
           << tags[1] << endln;
    return nullptr;
  }

  int argn = 4;
  const char *dmgtag = 0;
  double Ccd = 0.5;
  while (argn < argc) {
    const char *param = argv[argn];

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-dmg") == 0)    ||
        (strcmp(param, "-DMG") == 0))   {
      dmgtag = argv[++argn];
    } else if ((strcmp(param, "-couple") == 0) || 
               (strcmp(param, "-ccd") == 0)    ||
               (strcmp(param, "-Ccd") == 0))   {
      Ccd = std::stod(argv[++argn]);
      // opserr << "WARNING invalid baseTag uniaxialMaterial " WRAPPER_CMD ;
    } else {
      opserr << "WARNING invalid option: " << param
             << " in uniaxialMaterial '" WRAPPER_CMD "' with tag: '" 
             << tags[0] << "'"
             << endln;
      return nullptr;
    }
    argn++;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial,
                                              minStrain, maxStrain);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type " WRAPPER_CMD << endln;
    return nullptr;
  }
  theMaterial->setCoupling(Ccd);

  if (dmgtag){
    if (theMaterial->setDamageWrapper(G3_getInterpreter(rt), dmgtag) > 0)
      opserr << "#Set damage wrapper '" << dmgtag << "'\n";
  }

  // return G3_addUniaxialMaterial(rt, theMaterial);
  return theMaterial;
}

