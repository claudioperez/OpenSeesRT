#include <OPS_Stream.h>

void OPS_printCrdTransf(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"crdTransformations\": [\n";        
    MapOfTaggedObjectsIter theObjects = theCrdTransfObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theCrdTransfObjects.getNumComponents();    
    while ((theObject = theObjects()) != 0) {
      CrdTransf *theTransf = (CrdTransf *)theObject;
      theTransf->Print(s, flag);
      if (count < numComponents-1)
        s << ",\n";
      count++;      
    }
    s << "\n\t\t]";
  }
}

void OPS_printNDMaterial(OPS_Stream &s, int flag) {
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\"ndMaterials\": [\n";
        MapOfTaggedObjectsIter theObjects = theNDMaterialObjects.getIter();
        theObjects.reset();
        TaggedObject *theObject;
        int count = 0;
        int numComponents = theNDMaterialObjects.getNumComponents();
        while ((theObject = theObjects()) != 0) {
            NDMaterial *theMaterial = (NDMaterial *)theObject;
            theMaterial->Print(s, flag);
            if (count < numComponents - 1)
                s << ",\n";
            count++;
        }
        s << "\n\t\t]";
    }
}

void OPS_printSectionForceDeformation(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"sections\": [\n";    
    MapOfTaggedObjectsIter theObjects = theSectionForceDeformationObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theSectionForceDeformationObjects.getNumComponents();
    while ((theObject = theObjects()) != 0) {
      SectionForceDeformation *theSection = (SectionForceDeformation *)theObject;
      theSection->Print(s, flag);
      if (count < numComponents-1)
              s << ",\n";
      count++;
    }
    s << "\n\t\t]";
  }
}

void OPS_printUniaxialMaterial(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"uniaxialMaterials\": [\n";        
    MapOfTaggedObjectsIter theObjects = theUniaxialMaterialObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theUniaxialMaterialObjects.getNumComponents();    
    while ((theObject = theObjects()) != 0) {
      UniaxialMaterial *theMaterial = (UniaxialMaterial *)theObject;
      theMaterial->Print(s, flag);
      if (count < numComponents-1)
        s << ",\n";
      count++;      
    }
    s << "\n\t\t]";
  }
}

// Talledo Start
int
printModelGID(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char **argv)
{
  // This function print's a file with node and elements in a format useful for
  // GID
  int res = 0;
  bool hasLinear = 0;
  bool hasTri3 = 0;
  bool hasQuad4 = 0;
  bool hasQuad8 = 0;
  bool hasQuad9 = 0;
  bool hasBrick = 0;
  int startEle = 1;
  int endEle = 1;
  int eleRange = 0;
  int i = 2;

  FileStream outputFile;
  OPS_Stream *output = &opserr;

  if (argc < 2) {
    opserr << "WARNING printGID fileName? - no filename supplied\n";
    return TCL_ERROR;
  }
  openMode mode = OVERWRITE;
  if (argc >= 3) {
    if (strcmp(argv[i], "-append") == 0) {
      mode = APPEND;
      i++;
    }
    if (strcmp(argv[i], "-eleRange") == 0) {
      // opserr<<"WARNING:commands: eleRange defined"<<endln;
      eleRange = 1;
      if (Tcl_GetInt(interp, argv[i + 1], &startEle) != TCL_OK) {
        opserr << "WARNING print node failed to get integer: " << argv[i + 1]
               << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i + 2], &endEle) != TCL_OK) {
        opserr << "WARNING print node failed to get integer: " << argv[i + 2]
               << endln;
        return TCL_ERROR;
      }
      // opserr<<"startEle = "<<startEle<<" endEle = "<<endEle<<endln;
    }
  }

  if (outputFile.setFile(argv[1], mode) < 0) {
    opserr << "WARNING printGID " << argv[1] << " failed to set the file\n";
    return TCL_ERROR;
  }

  // Cycle over Elements to understand what type of elements are there
  ElementIter &theElements = theDomain.getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) {

    // Check type of Element with Number of Nodes
    // if 2 Nodes print the Element
    int nNode = theElement->getNumExternalNodes();
    if (nNode == 2) {
      hasLinear = 1;
    } else if (nNode == 4) {
      hasQuad4 = 1;
    } else if (nNode == 3) {
      hasTri3 = 1;
    } else if (nNode == 9) {
      hasQuad9 = 1;
    } else if (nNode == 8) {
      const char *name = theElement->getClassType();
      if (strcmp(name, "Brick") == 0) {
        hasBrick = 1;
      } else {
        hasQuad8 = 1;
      }
    }
  }
  // **** Linear Elements - 2 Nodes
  if (hasLinear == 1) {
    // Print HEADER
    outputFile << "MESH \"2NMESH\" dimension 3 ElemType Linear Nnode 2"
               << endln;
    outputFile << "#color 0 0 255" << endln << endln;

    // Print node coordinates
    outputFile << "Coordinates" << endln;
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << endln;
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << endln;
    }
    outputFile << "End coordinates" << endln << endln;

    // Print elements connectivity
    outputFile << "Elements" << endln;
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {
        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 2) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; i++) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; i++) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << endln;
        }
      }
    }
    outputFile << "End elements" << endln;
  }
  // **** Quadrilateral Elements - 4 Nodes
  if (hasQuad4 == 1) {
    // Print HEADER
    outputFile << "MESH \"4NMESH\" dimension 3 ElemType Quadrilateral Nnode 4"
               << endln;
    outputFile << "#color 0 255 0" << endln << endln;

    // Print node coordinates
    outputFile << "Coordinates" << endln;
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << endln;
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << endln;
    }
    outputFile << "End coordinates" << endln << endln;

    // Print elements connectivity
    outputFile << "Elements" << endln;
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 4) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; i++) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; i++) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << endln;
        }
      }
    }
    outputFile << "End elements" << endln;
  }
  // **** Triangular Elements - 3 Nodes
  if (hasTri3 == 1) {
    // Print HEADER
    outputFile << "MESH \"3NMESH\" dimension 3 ElemType Triangle Nnode 3"
               << endln;
    outputFile << "#color 0 255 0" << endln << endln;

    // Print node coordinates
    outputFile << "Coordinates" << endln;
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << endln;
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << endln;
    }
    outputFile << "End coordinates" << endln << endln;

    // Print elements connectivity
    outputFile << "Elements" << endln;
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 3 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 3) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; i++) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; i++) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << endln;
        }
      }
    }
    outputFile << "End elements" << endln;
  }
  // **** Quadrilateral Elements - 9 Nodes
  if (hasQuad9 == 1) {
    // Print HEADER
    outputFile << "MESH \"9NMESH\" dimension 3 ElemType Linear Nnode 9"
               << endln;
    outputFile << "#color 0 255 0" << endln << endln;

    // Print node coordinates
    outputFile << "Coordinates" << endln;
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << endln;
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << endln;
    }
    outputFile << "End coordinates" << endln << endln;

    // Print elements connectivity
    outputFile << "Elements" << endln;
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 9) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; i++) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; i++) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << endln;
        }
      }
    }
    outputFile << "End elements" << endln;
  }
  // **** Hexahedra Elements - 8 Nodes
  if (hasBrick == 1) {
    // Print HEADER
    outputFile << "MESH \"8NMESH\" dimension 3 ElemType Hexahedra Nnode 8"
               << endln;
    outputFile << "#color 255 0 0" << endln << endln;

    // Print node coordinates
    outputFile << "Coordinates" << endln;
    NodeIter &theNodes = theDomain.getNodes();
    MeshRegion *myRegion = theDomain.getRegion(0);
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << endln;
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << endln;
    }
    outputFile << "End coordinates" << endln << endln;

    // Print elements connectivity
    outputFile << "Elements" << endln;
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 8) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; i++) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; i++) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << endln;
        }
      }
    }
    outputFile << "End elements" << endln;
  }

  outputFile.close();
  return res;
}
// Talledo End
