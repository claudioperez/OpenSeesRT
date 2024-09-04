
int
TclBasicBuilder_addFourNodeQuad(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;


  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? "
              "thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int FourNodeQuadId, iNode, jNode, kNode, lNode, matID;
  double thickness = 1.0;
  double p = 0.0;   // uniform normal traction (pressure)
  double rho = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadId) != TCL_OK) {
    opserr << "WARNING invalid FourNodeQuad eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << "WARNING invalid kNode\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << "WARNING invalid lNode\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
    opserr << "WARNING invalid thickness\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  TCL_Char *type = argv[6 + argStart];

  if (Tcl_GetInt(interp, argv[7 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) > 11) {
    if (Tcl_GetDouble(interp, argv[8 + argStart], &p) != TCL_OK) {
      opserr << "WARNING invalid pressure\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9 + argStart], &rho) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;

  // now create the FourNodeQuad and add it to the Domain
  FourNodeQuad *theFourNodeQuad =
      new FourNodeQuad(FourNodeQuadId, iNode, jNode, kNode, lNode, *theMaterial,
                       type, thickness, p, rho, b1, b2);

  if (builder->getDomain()->addElement(theFourNodeQuad) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "FourNodeQuad element: " << FourNodeQuadId << "\n";
    delete theFourNodeQuad;
    return TCL_ERROR;
  }

  return TCL_OK;
}
