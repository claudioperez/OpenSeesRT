#include <g3_api.h>
#include <stdio.h>
#include <g3_api.h>
#include <DegradingUniaxialWrapper.hh>

#define WRAPPER_CMD "FedeasDegradingWrapper"

// TODO: change to TclSafeBuildObj_
int
TclSafeBuilder_addFedeasWrapper(ClientData clientData, G3_Runtime *rt,
                                  int argc, TCL_Char **argv)
{
  // Pointer to a uniaxial material that will be returned
  DegradingUniaxialWrapper *theMaterial = 0;
  UniaxialMaterial *theWrappedMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain =  1.0e16;
  int tags[2];

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-min $minStrain> <-max $maxStrain>"
           << endln;
    return 0;
  }

  // Get wrapper tag
  if (Tcl_GetInt(rt, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return 0;
  }
  // Get base tag
  if (Tcl_GetInt(rt, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return 0;
  }

  // Get base material
  theWrappedMaterial = G3_getUniaxialMaterialInstance(rt, tags[1]);
  if (theWrappedMaterial == 0) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag" WRAPPER_CMD " tag: "
           << tags[1] << endln;
    return 0;
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
      return 0;
    }
    argn++;
  }

// Parsing was successful, allocate the material
  theMaterial = new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial,
                                              minStrain, maxStrain);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type " WRAPPER_CMD << endln;
    return 0;
  }
  theMaterial->setCoupling(Ccd);
  if (dmgtag){
    if (theMaterial->setDamageWrapper(rt, dmgtag) > 0)
      opserr << "#Set damage wrapper '" << dmgtag << "'\n";
  }

  return G3_addUniaxialMaterial(rt, theMaterial);
  // return theMaterial;
}

