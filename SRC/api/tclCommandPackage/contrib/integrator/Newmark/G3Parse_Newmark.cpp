#include <tcl.h>
#include <g3_api.h>
#include <InputAPI.h>
#include <G3_Logging.h>
// #include <G3_Runtime.h>
// #include <SRC/analysis/integrator/Newmark.h>
#include <Newmark.h>

TransientIntegrator*
G3Parse_newNewmarkIntegrator(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** argv)
{
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = nullptr;

  if (argc != 4 && argc != 6) {
    opserr << G3_ERROR_PROMPT << " incorrect number of args want Newmark $gamma $beta "
              "<-form $typeUnknown>\n";
    opserr << "        got ";
    for (int i=0; i<argc; i++) opserr << argv[i] << ",";
    opserr << "\n";
    return 0;
  }

  int argi = 2;
  int dispFlag = 1;
  double dData[2];
  int numData = 2;
  for (; argi < 4; argi++)
    if (Tcl_GetDouble(interp, argv[argi], &dData[argi-2]) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid arg at position '" << argi << "'. Expected:\n";
      opserr << "\tintegrator Newmark $gamma $beta <-form $typeUnknown>\n";
      opserr << "  but got '" << argv[argi] << "'.\n";
      return nullptr;
    }
  // opserr << "Newmark(" << dData[0] << ", " << dData[1] << ")\n";

  if (argc == 4)
    return new Newmark(dData[0], dData[1]);

  else {
    const char *nextString = argv[argi];
    if (strcmp(nextString, "-form") == 0) {
      nextString = argv[++argi];
      if ((nextString[0] == 'D') || (nextString[0] == 'd'))
        dispFlag = 1;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        dispFlag = 3;
      else if ((nextString[0] == 'V') || (nextString[0] == 'v'))
        dispFlag = 2;
    }
    return new Newmark(dData[0], dData[1], dispFlag);
  }

  return nullptr;
}
