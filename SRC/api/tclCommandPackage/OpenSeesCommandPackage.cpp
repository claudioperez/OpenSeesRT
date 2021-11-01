/* Claudio Perez */
#include <tcl.h>

int myCommands(Tcl_Interp *interp);
int OpenSeesAppInit(Tcl_Interp *interp);

// Error streams
#include <handler/OPS_Stream.h>
#include <StandardStream.h>
// Create global error stream
/*
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
void *simulationInfo=0;
*/

extern "C" {
static int
Hello_Cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  Tcl_SetObjResult(interp, Tcl_NewStringObj("Hello, World!", -1));
  return TCL_OK;
}

/*
 * Hello_Init -- Called when Tcl loads your extension.
 */
int DLLEXPORT
Openseescommandpackage_Init(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }

  if (Tcl_PkgProvide(interp, "OpenSeesCommandPackage", "0.0.1") == TCL_ERROR) {
    return TCL_ERROR;
  }

  Tcl_CreateObjCommand(interp, "hello", Hello_Cmd, NULL, NULL);

  OpenSeesAppInit(interp);
  myCommands(interp);
  return TCL_OK;
}
} // extern "C"

