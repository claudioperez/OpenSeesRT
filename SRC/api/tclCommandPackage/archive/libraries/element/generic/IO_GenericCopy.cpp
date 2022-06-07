

#include <g3_api.h>


#include <SRC/element/generic/GenericCopy.h>
void *OPS_GenericCopy()
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element genericCopy eleTag -node Ndi ... -src srcTag\n";
    return 0;
  }

  // tags
  int tag;
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING: invalid tag\n";
    return 0;
  }

  // nodes
  const char *type = OPS_GetString();
  if (strcmp(type, "-node") != 0) {
    opserr << "WARNING expecting -node Ndi Ndj ...\n";
    return 0;
  }
  ID nodes(32);
  int numNodes = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    int node;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &node) < 0) {
      break;
    }
    nodes(numNodes++) = node;
  }
  nodes.resize(numNodes);

  // source element
  int srcTag;
  numdata = 1;
  type = OPS_GetString();
  if (strcmp(type, "-src") != 0) {
    opserr << "WARNING expecting -src srcTag\n";
    return 0;
  }
  if (OPS_GetIntInput(&numdata, &srcTag) < 0) {
    opserr << "WARNING: invalid srcTag\n";
    return 0;
  }

  // create object
  Element *theEle = new GenericCopy(tag, nodes, srcTag);

  return theEle;
}
