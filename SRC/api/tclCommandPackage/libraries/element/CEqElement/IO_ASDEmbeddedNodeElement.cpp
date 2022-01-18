

#include <g3_api.h>


#include <SRC/element/CEqElement/ASDEmbeddedNodeElement.h>
void *OPS_ASDEmbeddedNodeElement(void)
{
  static bool first_done = false;
  if (!first_done) {
    opserr << "Using ASDEmbeddedNodeElement - Developed by: Massimo Petracca, "
              "Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  const char *descr = "Want: element ASDEmbeddedNodeElement $tag $Cnode "
                      "$Rnode1 $Rnode2 $Rnode3 <$Rnode4> <-rot> <-K $K>\n";

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 5) {
    opserr << "ASDEmbeddedNodeElement ERROR : Few arguments:\n" << descr;
    return 0;
  }

  // mandatory parameters
  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "ASDEmbeddedNodeElement ERROR: Invalid integer mandatory values: "
              "element ASDEmbeddedNodeElement wants at least 5 integer "
              "parameters\n"
           << descr;
    return 0;
  }

  // parse optional parameters
  bool rot = false;
  int N4 = 0;
  bool has_N4 = false;
  double K = 1.0e18;
  for (int i = 5; i < numArgs; i++) {
    const char *what = OPS_GetString();
    if (strcmp(what, "-rot") == 0) {
      rot = true;
    } else if (strcmp(what, "-K") == 0) {
      if (i == numArgs - 1) {
        opserr << "ASDEmbeddedNodeElement ERROR: The -K keyword should be "
                  "followed by a floating point number.\n"
               << descr;
        return 0;
      }
      ++i;
      numData = 1;
      if (OPS_GetDouble(&numData, &K) != 0) {
        opserr << "ASDEmbeddedNodeElement ERROR invalid floating point number "
                  "for -K keyword.\n";
        return 0;
      }
    } else {
      // should be an integer for the only if i == 6
      if (i == 5) {
        try {
          N4 = std::stoi(std::string(what));
          has_N4 = true;
        } catch (...) {
          N4 = -1;
          has_N4 = false;
        }
      }
    }
  }

  // done
  if (has_N4)
    return new ASDEmbeddedNodeElement(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], N4, rot, K);
  else
    return new ASDEmbeddedNodeElement(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], rot, K);
}
