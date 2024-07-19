//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addTruss()
// command. 
//
#include <assert.h>
#include <string.h>
// OpenSees
#include <tcl.h>
#include <Domain.h>
#include <Logging.h>
#include <Parsing.h>
#include <FrameSection.h>
#include <UniaxialMaterial.h>
#include <BasicModelBuilder.h>
// Elements
#include <Truss.h>
#include <TrussSection.h>
#include <CorotTruss.h>
#include <CorotTrussSection.h>

//
//  element Truss        $tag $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag> <-useInitialDisp $flag>
//  element Truss        $tag $iNode $jNode $sectTag   <-rho $rho> <-cMass $flag> <-doRayleigh $flag>
//  element TrussSection $tag $iNode $jNode $sectTag   <-rho $rho> <-cMass $flag> <-doRayleigh $flag>
//
int
TclCommand_addTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
                         TCL_Char **argv)
{

  constexpr int eleArgStart = 1;

  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element truss eleTag? iNode? jNode? A? matTag?\n";
    opserr << " -or- element truss eleTag? iNode? jNode? sectTag?" << "\n";
    return TCL_ERROR;
  }    

  int ndm = builder->getNDM();

  // get the id and end nodes 
  int trussId, iNode, jNode, matID;
  double rho = 0.0;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &trussId) != TCL_OK) {
    opserr << "WARNING invalid truss eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "truss element: " << trussId << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "truss element: " << trussId << "\n";
     return TCL_ERROR;
  }

  for (int i = 4+eleArgStart; i < argc; i++) {
    if (i+1 < argc && strcmp(argv[i], "-rho") == 0) {
      if (Tcl_GetDouble(interp, argv[i+1], &rho) != TCL_OK) {
        opserr << "WARNING invalid rho\n";
        opserr << "truss element: " << trussId << "\n";
        return TCL_ERROR;
      }
      argc -= 2;
      break;
    }
  }

  // UniaxialMaterial form
  if ((argc-eleArgStart) == 6) {
      double A;
      if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
          opserr << "WARNING invalid A\n";
          opserr << "truss element: " << trussId << "\n";
          return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[5+eleArgStart], &matID) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "truss element: " << trussId << "\n";
          return TCL_ERROR;
      }
      
      UniaxialMaterial *theMaterial = builder->getTypedObject<UniaxialMaterial>(matID); 
      if (theMaterial == 0)
          return TCL_ERROR;

      // now create the truss and add it to the Domain
      Element *theTruss = nullptr;
      if (strcmp(argv[eleArgStart],"corotTruss") == 0)
          theTruss = new CorotTruss(trussId,ndm,iNode,jNode,*theMaterial,A,rho);
      else
          theTruss = new Truss(trussId,ndm,iNode,jNode,*theMaterial,A,rho);


      if (builder->getDomain()->addElement(theTruss) == false) {
          opserr << "WARNING could not add element to the domain\n";
          opserr << "truss element: " << trussId << "\n";
          delete theTruss;
          return TCL_ERROR;
      }

  // TrussSection
  } else {
      int sectID;
      if (Tcl_GetInt(interp, argv[4+eleArgStart], &sectID) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "truss element: " << trussId << "\n";
          return TCL_ERROR;
      }
      
      SectionForceDeformation *theSection = builder->getTypedObject<FrameSection>(sectID);

      if (theSection == nullptr)
          return TCL_ERROR;
  
      // now create the truss and add it to the Domain
      Element *theTruss = 0;
      if (strcmp(argv[eleArgStart],"corotTruss") == 0)
          theTruss = new CorotTrussSection(trussId,ndm,iNode,jNode,*theSection,rho);
      else
          theTruss = new TrussSection(trussId,ndm,iNode,jNode,*theSection,rho);

      if (builder->getDomain()->addElement(theTruss) == false) {
          opserr << "WARNING could not add element to the domain\n";
          opserr << "truss element: " << trussId << "\n";
          delete theTruss;
          return TCL_ERROR;
      }
  }
  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



