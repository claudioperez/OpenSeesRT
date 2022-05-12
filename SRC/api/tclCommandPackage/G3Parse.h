#include <G3_Runtime.h>
#include <tcl.h>

#define G3_Char TCL_Char

// #define OPS_GetIntInput(ndat, dat)  (rt->interp->OPS_GetIntInput((ndat), (dat)))
// #define OPS_GetDoubleInput (rt->getDoubleInput)

#define G3Parse_getInt(bld, arg, adr)    Tcl_GetInt((bld)->m_interp, (arg), (adr))
#define G3Parse_getDouble(bld, arg, adr) Tcl_GetDouble((bld)->m_interp, (arg), (adr))
// #define G3Parse_AppendResult(interp, 


