#ifndef OPS_RUNTIME_API_H
#define OPS_RUNTIME_API_H

class G3_Runtime;
#ifndef OPS_USE_RUNTIME
#  define OPS_ADD_RUNTIME_VPV(func) (func)()
#  define OPS_CALL_RUNTIME_VPV(func) (func)()
#  define OPS_ADD_RUNTIME_IXV(func) (func)()
#  define OPS_ADD_RUNTIME_VXV(func) (func)()
#  define OPS_DECL_RUNTIME_VPID(func, arg) func(arg)
#  define OPS_DECL_RUNTIME(func, ...) (func)(__VA_ARGS__)

#else
#  define OPS_ADD_RUNTIME_VPV(func) (func)([[maybe_unused]] G3_Runtime *rt, int _argc, const char **const argv)
#  define OPS_CALL_RUNTIME_VPV(func) (func)(rt, 0, 0)
#  define OPS_ADD_RUNTIME_IXV(func) (func)([[maybe_unused]] G3_Runtime *rt)
#  define OPS_ADD_RUNTIME_VXV(func) (func)([[maybe_unused]] G3_Runtime *rt)
#  define OPS_DECL_RUNTIME_VPID(func, arg_t) (func)([[maybe_unused]] G3_Runtime *rt, arg_t)
#  define OPS_DECL_RUNTIME(func, ...) (func)([[maybe_unused]] G3_Runtime *rt, __VA_ARGS__)
#  undef OPS_GetDomain
#  include "internal/g3_api.h"
   typedef void *OPS_Routine(G3_Runtime* , int, const char** const);
#endif

#endif // OPS_RUNTIME_API_H
