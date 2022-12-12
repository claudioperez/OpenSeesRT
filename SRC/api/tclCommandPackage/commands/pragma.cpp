/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <tcl.h>
#include <string.h>


int
TclObjCommand_pragma([[maybe_unused]] ClientData clientData, 
                     Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  if (objc == 1)
    return TCL_OK;

  if (objc == 2)
    return TCL_OK;
  
  int argi = 1;
  if (strcmp(Tcl_GetString(objv[argi++]), "analysis") == 0) {
    if (strcmp(Tcl_GetString(objv[argi]), "off") == 0) {
      Tcl_Eval(interp,
        "proc loadConst {args} {}\n"
        "proc wipeAnalysis	{args} {}\n"
        "proc constraints {args} {}\n"
        "proc numberer {args} {}\n"
        "proc system {args} {}\n"
        "proc test {args} {}\n"
        "proc algorithm {args} {}\n"
        "proc integrator {args} {}\n"
        "proc analysis {args} {}\n"
        "proc analyze {args} {}\n"
        "namespace eval opensees::pragma {set analysis off}\n"
      );
      return TCL_OK;
    } 
  }
  return TCL_OK;
}

