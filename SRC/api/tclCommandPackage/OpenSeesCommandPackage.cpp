/* Claudio Perez */
#include <g3_api.h>

int myCommands(G3_Runtime *rt);
int OpenSeesAppInit(G3_Runtime *rt);

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
Hello_Cmd(ClientData cdata, G3_Runtime *rt, int objc, Tcl_Obj *const objv[])
{
  Tcl_SetObjResult(rt, Tcl_NewStringObj("Hello, World!", -1));
  return TCL_OK;
}

/*
 * Hello_Init -- Called when Tcl loads your extension.
 */
int DLLEXPORT
Openseescommandpackage_Init(G3_Runtime *rt)
{
  if (Tcl_InitStubs(rt, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }

  if (Tcl_PkgProvide(rt, "OpenSeesCommandPackage", "0.0.1") == TCL_ERROR) {
    return TCL_ERROR;
  }

  Tcl_CreateObjCommand(rt, "hello", Hello_Cmd, NULL, NULL);

  OpenSeesAppInit(rt);
  myCommands(rt);
  return TCL_OK;
}
} // extern "C"

