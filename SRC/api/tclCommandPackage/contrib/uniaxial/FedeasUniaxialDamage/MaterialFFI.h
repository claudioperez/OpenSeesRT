#pragma once
#include "ArgumentFFI.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ISW_COMMIT
#undef ISW_COMMIT
#endif
#ifdef ISW_DELETE
#undef ISW_DELETE
#endif

enum ISW_Task {
  ISW_COMMIT = 1<<0,
  ISW_DELETE = 1<<1,
  ISW_UPDATE = 1<<2,
  ISW_MALLOC = 1<<3,
  ISW_CREATE = 1<<4,

  ISW_RETURN_RESIDUAL = 1<<5,
  ISW_RETURN_TANGENT  = 1<<6,
  ISW_RETURN_OTHER    = 1<<7,
  
  ISW_UPDATE_InitialTangent=1<<8,
  ISW_UPDATE_CurrentTangent=1<<9
};

const int ISW_ACTION = ISW_COMMIT 
               | ISW_UPDATE 
               | ISW_MALLOC 
               | ISW_DELETE 
               | ISW_CREATE;

const int ISW_MODIFY = ISW_UPDATE_InitialTangent
               | ISW_UPDATE_CurrentTangent;


struct Runtime {
  int ndm, ndf;
  double time, time_step;
  G3_Arg *const objv[];

//   int    argc;
//   char** argv;
};

struct MaterialRoutine;

typedef int (MaterialFunction)(
       enum           ISW_Task,
       void*          runtime,
       struct MaterialRoutine*,
       struct MaterialResponse*);

struct MaterialRoutine {
  MaterialFunction *routine;

  double *variable,
         *velocity,
         *acceleration;

  // Parsing
  const char  **argv;
  int  argc;

  void   *data, *instance;
  size_t data_size;
};

struct MaterialResponse {
  double *response,
         *residual,
         *jacobian;
  size_t response_size;
};

#ifdef __cplusplus
}
#endif
