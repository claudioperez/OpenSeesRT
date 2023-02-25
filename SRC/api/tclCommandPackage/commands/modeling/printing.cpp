/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#include <assert.h>
#include <g3_api.h>
#include <FileStream.h>
#include <SimulationInformation.h>

#include <BasicModelBuilder.h>
#include <Storage/G3_TableIterator.h>

#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>

#include <Pressure_Constraint.h>
#include <Element.h>
#include <ElementIter.h>

#include <Node.h>
#include <NodeIter.h>

extern  SimulationInformation simulationInfo;

int printElement(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char **argv, OPS_Stream &output);

int printNode(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv, OPS_Stream &output);

int printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char **argv, OPS_Stream &output);

int printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char **argv, OPS_Stream &output);


static int
printRegistry(BasicModelBuilder* builder, TCL_Char* type, int flag, OPS_Stream *output)
{
  G3_TableIterator iter = builder->iterate(type);

  while (G3_NextTableEntry(&iter)) {
    ((TaggedObject*)iter.value)->Print(*output, flag);
    (*output) << ",\n";
  }
  
  return TCL_OK;
}

static void
printDomain(OPS_Stream &s, BasicModelBuilder* builder, int flag) 
{

  Domain* theDomain = builder->getDomain();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {

    s << "\t\"properties\": {\n";

#if 0
    OPS_printNDMaterial(s, flag);
    s << ",\n";
    OPS_printSectionForceDeformation(s, flag);
    s << ",\n";   
    OPS_printCrdTransf(s, flag);      
    s << ",\n";
#endif
    //
    s << "\t\t\"uniaxialMaterials\": [\n";        
    printRegistry(builder, "UniaxialMaterial", flag, &s);
    s << "\n\t\t]";
    s << ",\n";
    //
    s << "\t\t\"crdTransformations\": [\n";        
    printRegistry(builder, "CoordinateTransform", flag, &s);
    s << "\n\t\t]";
    //

    s << "\n\t},\n";
    s << "\t\"geometry\": {\n";

    int numToPrint = theDomain->getNumNodes();
    NodeIter &theNodess = theDomain->getNodes();
    Node *theNode;
    int numPrinted = 0;
    s << "\t\t\"nodes\": [\n";
    while ((theNode = theNodess()) != 0) {    
      theNode->Print(s, flag);
      numPrinted += 1;
      if (numPrinted < numToPrint)
        s << ",\n";
      else
        s << "\n\t\t],\n";
    }


    Element *theEle;
    ElementIter &theElementss = theDomain->getElements();
    numToPrint = theDomain->getNumElements();
    numPrinted = 0;
    s << "\t\t\"elements\": [\n";
    while ((theEle = theElementss()) != 0) {
      theEle->Print(s, flag);
      numPrinted += 1;
      if (numPrinted < numToPrint)
        s << ",\n";
      else
        s << "\n\t\t]\n";
      }

      s << "\t}\n";
      s << "}\n";
      s << "}\n";

      return;
  }
      
#if 0 
  s << "Current Domain Information\n";
  s << "\tCurrent Time: " << theDomain->getCurrentTime();
  // s << "\ntCommitted Time: " << committedTime << endln;    
  s << "NODE DATA: NumNodes: " << theDomain->getNumNodes() << "\n";  
  theNodes->Print(s, flag);
  
  s << "ELEMENT DATA: NumEle: " << theElements->getNumComponents() << "\n";
  theElements->Print(s, flag);
  
  s << "\nSP_Constraints: numConstraints: " << theSPs->getNumComponents() << "\n";
  theSPs->Print(s, flag);
  
  s << "\nPressure_Constraints: numConstraints: " << thePCs->getNumComponents() << "\n";
  thePCs->Print(s, flag);
  
  s << "\nMP_Constraints: numConstraints: " << theMPs->getNumComponents() << "\n";
  theMPs->Print(s, flag);
  
  s << "\nLOAD PATTERNS: numPatterns: " << theLoadPatterns->getNumComponents() << "\n\n";
  theLoadPatterns->Print(s, flag);
  
  s << "\nPARAMETERS: numParameters: " << theParameters->getNumComponents() << "\n\n";
  theParameters->Print(s, flag);
#endif
}

int
printModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  assert(clientData != nullptr);

  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain * domain = builder->getDomain();

  int currentArg = 1;
  int res = 0;

  int flag = OPS_PRINT_CURRENTSTATE;

  FileStream outputFile;
  OPS_Stream *output = &opserr;
  bool done = false;

  // if just 'print' then print out the entire domain
  if (argc == currentArg) {
    opserr << domain;
    return TCL_OK;
  }

  while (done == false) {
    // if 'print ele i j k..' print out some elements
    if ((strcmp(argv[currentArg], "-ele") == 0) ||
        (strcmp(argv[currentArg], "ele") == 0)) {
      currentArg++;
      res = printElement((ClientData)domain, interp, argc - currentArg, argv + currentArg, *output);
      done = true;
    }
    // if 'print node i j k ..' print out some nodes
    else if ((strcmp(argv[currentArg], "-node") == 0) ||
             (strcmp(argv[currentArg], "node") == 0)) {
      currentArg++;
      res = printNode((ClientData)domain, interp, argc - currentArg, argv + currentArg,
                      *output);
      done = true;
    }

    else if ((strcmp(argv[currentArg], "-registry") == 0)) {
      currentArg++;
      res = printRegistry((BasicModelBuilder*)clientData, argv[currentArg++], flag, output);
      done = true;
    }

    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[currentArg], "integrator") == 0) ||
             (strcmp(argv[currentArg], "-integrator") == 0)) {
      currentArg++;
      res = printIntegrator((ClientData)domain, interp, argc - currentArg,
                            argv + currentArg, *output);
      done = true;
    }

    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[currentArg], "algorithm") == 0) ||
             (strcmp(argv[currentArg], "-algorithm") == 0)) {
      currentArg++;
      res = printAlgorithm((ClientData)domain, interp, argc - currentArg,
                           argv + currentArg, *output);
      done = true;
    }
    else if ((strcmp(argv[currentArg], "-JSON") == 0)) {
      currentArg++;
      flag = OPS_PRINT_PRINTMODEL_JSON;
    }

    else {
      if ((strcmp(argv[currentArg], "file") == 0) ||
          (strcmp(argv[currentArg], "-file") == 0))
        currentArg++;

      openMode mode = APPEND;
      if (flag == OPS_PRINT_PRINTMODEL_JSON)
        mode = OVERWRITE;
      if (outputFile.setFile(argv[currentArg], mode) != 0) {
        opserr << "print <filename> .. - failed to open file: "
               << argv[currentArg] << endln;
        return TCL_ERROR;
      }
      currentArg++;

      // if just 'print <filename>' then print out the entire domain to eof
      if (argc == currentArg) {
        if (flag == OPS_PRINT_PRINTMODEL_JSON)
          simulationInfo.Print(outputFile, flag);

        // domain->Print(outputFile, flag);
        printDomain(outputFile, builder, flag);
        return TCL_OK;
      }

      output = &outputFile;
    }
  }

  // close the output file
  outputFile.close();
  return res;
}

int
printElement(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char **argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;

  int flag   = 0; // default flag sent to a nodes Print() method
  int eleArg = 0;

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) ||
      (strcmp(argv[0], "-flag")) == 0) { // get the specified flag
    if (argc < 2) {
      opserr << "WARNING print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print ele failed to get integer flag: \n";
      opserr << argv[eleArg] << endln;
      return TCL_ERROR;
    }
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default
  if (argc == eleArg) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output, flag);
    return TCL_OK;

  } else {

    // otherwise print out the specified nodes i j k .. with flag
    int numEle = argc - eleArg;
    ID *theEle = new ID(numEle);
    for (int i = 0; i < numEle; i++) {
      int eleTag;
      if (Tcl_GetInt(interp, argv[i + eleArg], &eleTag) != TCL_OK) {
        opserr << "WARNING print -ele failed to get integer: " << argv[i]
               << endln;
        return TCL_ERROR;
      }
      (*theEle)(i) = eleTag;
    }

    the_domain->Print(output, 0, theEle, flag);
    delete theEle;
  }

  return TCL_OK;
}

// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// Parameters
//   nodeArg: integer equal to arg count to node plus 1
//   output:  output stream to which the results are sent
//
int
printNode(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv,
          OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method
  int nodeArg = 0;

  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != nullptr)
      theNode->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) || (strcmp(argv[0], "-flag") == 0)) {
    // get the specified flag
    if (argc <= nodeArg) {
      opserr << "WARNING print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print node failed to get integer flag: \n";
      opserr << argv[nodeArg] << endln;
      return TCL_ERROR;
    }
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag'
  //     print out all the nodes in the domain with flag
  if (nodeArg == argc) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != nullptr)
      theNode->Print(output, flag);
    return TCL_OK;
  } else {
    // otherwise print out the specified nodes i j k .. with flag
    int numNodes = argc - nodeArg;
    ID *theNodes = new ID(numNodes);
    for (int i = 0; i < numNodes; i++) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[nodeArg], &nodeTag) != TCL_OK) {
        opserr << "WARNING print node failed to get integer: " << argv[nodeArg]
               << endln;
        return TCL_ERROR;
      }
      (*theNodes)(i) = nodeTag;
      nodeArg++;
    }

    domain->Print(output, theNodes, 0, flag);
    delete theNodes;
  }

  return TCL_OK;
}

