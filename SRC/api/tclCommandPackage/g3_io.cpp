
int
G3_Interpreter::OPS_GetIntInput(int *numData, int *data)
{
  int size = *numData;

  for (int i = 0; i < size; i++) {
    if ((currentArg >= maxArg) ||
        (Tcl_GetInt(theInterp, currentArgv[currentArg], &data[i]) != TCL_OK)) {
      return -1;
    } else
      currentArg++;
  }

  return 0;
}
