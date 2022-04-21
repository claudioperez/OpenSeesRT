/* Claudio Perez */
#include <g3_api.h>
#undef G3_Runtime
#include "G3_Runtime.h"

extern int myCommands(Tcl_Interp *interp);
extern int OpenSeesAppInit(Tcl_Interp *interp);
extern int init_g3_tcl_utils(Tcl_Interp*);

// Error streams
#include "streams/G3_Logging.h"
#include <handler/OPS_Stream.h>
#include <StandardStream.h>
#include <unistd.h>
// Create global error stream


extern "C" {

/*
 * Called when Tcl loads the extension.
 */
int DLLEXPORT
Openseesrt_Init(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }

  if (Tcl_PkgProvide(interp, "OpenSeesCommandPackage", "0.0.1") == TCL_ERROR) {
    return TCL_ERROR;
  }

  Tcl_SetAssocData(interp, "G3_Runtime", NULL, (ClientData)(new G3_Runtime{interp}));

  OpenSeesAppInit(interp);
  myCommands(interp);
  init_g3_tcl_utils(interp);
  if (isatty(STDERR_FILENO))
    G3_setStreamColor(nullptr, G3_Warn, 1);
  return TCL_OK;
}

} // extern "C"

