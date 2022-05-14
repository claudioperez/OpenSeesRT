#ifndef G3PARSE_H
#include <G3_Runtime.h>
#include <tcl.h>

#define G3_Char TCL_Char

// #define OPS_GetIntInput(ndat, dat)  (rt->interp->OPS_GetIntInput((ndat), (dat)))
// #define OPS_GetDoubleInput (rt->getDoubleInput)

#define G3Parse_getInt(bld, arg, adr)    Tcl_GetInt((bld)->m_interp, (arg), (adr))
#define G3Parse_getDouble(bld, arg, adr) Tcl_GetDouble((bld)->m_interp, (arg), (adr))
// #define G3Parse_AppendResult(interp, 


static void printCommand(int argc, TCL_Char **argv) {
  opserr << "Input command: ";
  for (int i = 0; i < argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
}
#endif // G3PARSE_H
