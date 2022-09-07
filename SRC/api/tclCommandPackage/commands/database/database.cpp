
extern int TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv, Domain &theDomain,
                          FEM_ObjectBroker &theBroker);

int
addDatabase(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char **argv)
{
  return TclAddDatabase(clientData, interp, argc, argv, theDomain, theBroker);
}

