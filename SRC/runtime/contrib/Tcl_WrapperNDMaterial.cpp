
NDMaterial *
Tcl_addWrapperNDMaterial(matObj *theMat, ClientData clientData,
                         Tcl_Interp *interp, int argc, TCL_Char ** const argv,
                         TclBasicBuilder *theModelbuilder)
{
  theInterp = interp;

  // theModelBuilder = builder;
  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

  // get the current load factor
  static modelState theModelState;
  if (theDomain != 0) {
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;
    theModelState.time = time;
    theModelState.dt = dt;
  }

  // invoke the mat function with isw = 0
  int isw = ISW_INIT;
  int result = 0;
  theMat->matFunctPtr(theMat, &theModelState, 0, 0, 0, &isw, &result);
  int matType = theMat->matType; // GR added to support material

  if (result != 0 ||
      (matType != OPS_PLANESTRESS_TYPE && matType != OPS_PLANESTRAIN_TYPE &&
       matType != OPS_THREEDIMENSIONAL_TYPE)) {
    opserr << "Tcl_addWrapperNDMaterial - failed in element function " << result
           << endln;
    return 0;
  }

  WrapperNDMaterial *theMaterial =
      new WrapperNDMaterial(argv[1], theMat, theMat->matType);

  return theMaterial;
}

