
#include <tcl.h>
#include <string.h>

#include <unordered_map.h>

Tcl_CmdProc TclCommand_useUniaxialMaterial;

const struct {const char*name; const Tcl_CmdProc*func;} command_table[] = {
  {"UniaxialMaterial",  TclCommand_useUniaxialMaterial       },
};

TclCommand_invoke(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{

}


