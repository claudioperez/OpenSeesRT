
#define OPS_GetIntInput(ndat, dat)  (rt->interp->OPS_GetIntInput((ndat), (dat)))
#define OPS_GetDoubleInput (rt->getDoubleInput)

#define G3_GetInt(interp, arg, adr) Tcl_GetInt(rt->getInterp(), (arg), (adr))
#define G3_GetDouble(interp, arg, adr) Tcl_GetDouble(rt->getInterp(), (arg), (adr))
#define G3_AppendResult(interp, 


