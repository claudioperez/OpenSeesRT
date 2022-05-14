
#include <g3_api.h>
#include <G3Parse.h>
// #include <G3_Runtime.h>
// #include <SRC/analysis/integrator/Newmark.h>
#include <Newmark.h>

TransientIntegrator*
G3Parse_newNewmarkIntegrator(G3_Runtime*rt, int argc, G3_Char** argv)
{
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = nullptr;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want Newmark $gamma $beta "
              "<-form $typeUnknown>\n";
    opserr << "        got ";
    for (int i=0; i<argc; i++) opserr << argv[i] << ",";
    opserr << "\n";
    return 0;
  }

  int argi = 1;
  int dispFlag = 1;
  double dData[2];
  int numData = 2;
  for (; argi < 3; argi++)
    if (G3Parse_getDouble(rt, argv[argi], &dData[argi-1]) != TCL_OK) {
      opserr << "WARNING - invalid args want Newmark $gamma $beta <-form "
                "$typeUnknown>\n";
      return nullptr;
    }

  if (argc == 2)
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
