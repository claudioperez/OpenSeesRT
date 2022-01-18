

#include <g3_api.h>


#include <SRC/element/adapter/Adapter.h>
void *OPS_Adapter()
{
  int ndf = OPS_GetNDF();
  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element adapter eleTag -node Ndi Ndj ... -dof dofNdi -dof "
              "dofNdj ... -stif Kij ipPort <-ssl> <-udp> <-doRayleigh> <-mass "
              "Mij>\n";
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
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (OPS_GetIntInput(&numdata, &node) < 0) {
      if (numArgs > OPS_GetNumRemainingInputArgs()) {
        // move current arg back by one
        OPS_ResetCurrentInputArg(-1);
      }
      break;
    }
    nodes(numNodes++) = node;
  }
  nodes.resize(numNodes);

  // dofs
  int numDOF = 0;
  ID *dofs = new ID[numNodes];
  for (int i = 0; i < numNodes; i++) {
    type = OPS_GetString();
    if (strcmp(type, "-dof") != 0 && strcmp(type, "-dir") != 0) {
      opserr << "WARNING expecting -dof dofNd" << i + 1 << ", but got " << type
             << endln;
      return 0;
    }
    ID dofsi(ndf);
    int numDOFi = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      int dof;
      numdata = 1;
      int numArgs = OPS_GetNumRemainingInputArgs();
      if (OPS_GetIntInput(&numdata, &dof) < 0) {
        if (numArgs > OPS_GetNumRemainingInputArgs()) {
          // move current arg back by one
          OPS_ResetCurrentInputArg(-1);
        }
        break;
      }
      if (dof < 1 || ndf < dof) {
        opserr << "WARNING invalid dof ID\n";
        return 0;
      }
      dofsi(numDOFi++) = dof - 1;
      numDOF++;
    }
    dofsi.resize(numDOFi);
    dofs[i] = dofsi;
  }

  // stiffness matrix terms
  type = OPS_GetString();
  if (strcmp(type, "-stif") != 0 && strcmp(type, "-stiff") != 0) {
    opserr << "WARNING expecting -stif kij\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < numDOF * numDOF) {
    opserr << "WARNING wrong number of kij specified\n";
    return 0;
  }
  Matrix kb(numDOF, numDOF);
  numdata = 1;
  for (int i = 0; i < numDOF; i++) {
    for (int j = 0; j < numDOF; j++) {
      if (OPS_GetDoubleInput(&numdata, &kb(i, j)) < 0) {
        opserr << "WARNING invalid stiffness value\n";
        return 0;
      }
    }
  }
  // ipPort
  int ipPort;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &ipPort) < 0) {
    opserr << "WARNING: invalid ipPort\n";
    return 0;
  }

  // options
  int ssl = 0, udp = 0;
  int doRayleigh = 0;
  Matrix *mb = 0;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    return new Adapter(tag, nodes, dofs, kb, ipPort);
  }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    type = OPS_GetString();
    if (strcmp(type, "-ssl") == 0) {
      ssl = 1;
      udp = 0;
    } else if (strcmp(type, "-udp") == 0) {
      udp = 1;
      ssl = 0;
    } else if (strcmp(type, "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < numDOF * numDOF) {
        opserr << "WARNING wrong number of mij specified\n";
        return 0;
      }
      double mij;
      numdata = 1;
      mb = new Matrix(numDOF, numDOF);
      for (int i = 0; i < numDOF; i++) {
        for (int j = 0; j < numDOF; j++) {
          if (OPS_GetDoubleInput(&numdata, &mij) < 0) {
            opserr << "WARNING invalid damping value\n";
            delete mb;
            return 0;
          }
          (*mb)(i, j) = mij;
        }
      }
    }
  }

  // create object
  Element *theEle =
      new Adapter(tag, nodes, dofs, kb, ipPort, ssl, udp, doRayleigh, mb);

  // cleanup dynamic memory
  if (dofs != 0)
    delete[] dofs;
  if (mb != 0)
    delete mb;

  return theEle;
}
