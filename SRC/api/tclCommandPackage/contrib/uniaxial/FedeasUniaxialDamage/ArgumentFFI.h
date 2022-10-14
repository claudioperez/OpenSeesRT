#pragma once
#include <stdio.h>
#include <string.h>

#if 1
#  include <tcl.h>
#  define G3_ERROR TCL_ERROR
#  define G3_OK TCL_OK
   typedef Tcl_Obj G3_Arg;
#else
#  define G3_ERROR 0
#  define G3_OK 1
   typedef void* G3_Arg;
#endif
/*
int
G3_GetKey(void* interp, int argc, char** argv, char * key, int *loc)
{
  *loc = 0;
  int success = 1;
  while ((success = strcmp(key, argv[*loc++])) != 0)
    if (loc == argc)
      return G3_ERROR;
  return G3_OK;
}
*/

static int
G3_GetDouble([[maybe_unused]] void* interp, char *arg, double* addr) {
  if (sscanf(arg, "%lf", addr) == 1)
    return G3_OK;
  else
    return G3_ERROR;
}

static int
G3_GetInt([[maybe_unused]]void* interp, char *arg, int* addr) {
  if (sscanf(arg, "%d", addr) == 1)
    return G3_OK;
  else
    return G3_ERROR;
}

static void
G3_PrintCommand(int argc, char **argv)
{
  FILE* file = stderr;
  fprintf(file, "Command: ");
  for (int i = 0; i < argc; i++)
    fprintf(file, "%s ", argv[i]);
  fprintf(file, "\n");
}

