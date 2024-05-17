/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.15 $
// $Date: 2007-10-13 00:51:14 $
// $Source:
// /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/TclBeamWithHingesBuilder.cpp,v
// $

// File: ~/tcl/TclElmtBuilder.C
//
// Description: This file contains the implementation of the commands used
// to add sections and nonlinear frame elements to the model.
//
// Written: Remo M. de Souza
// Created: 08/99
// based on TclPlaneFrame.C by fmk and rms
//

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

#include <HingeMidpointBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>

#include <ElasticSection2d.h>
#include <ElasticSection3d.h>

#include <SectionForceDeformation.h>

#include <CrdTransf.h>

#include <runtime/BasicModelBuilder.h>

class TclBasicBuilder;
int
TclBasicBuilder_addBeamWithHinges(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv, Domain *theDomain,
                                  TclBasicBuilder *_builder)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  int NDM = builder->getNDM();
  int NDF = builder->getNDF();

  // Plane frame element
  if (NDM == 2 && NDF == 3) {
    if (argc < 13) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> "
                "<-iter maxIters tolerance>"
             << endln;
      return TCL_ERROR;
    }

    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, I;
    double massDens = 0.0;
    int numIters = 10;
    double tol = 1.0e-10;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &I) != TCL_OK) {
      opserr << "WARNING invalid I\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[12], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    bool isShear = false;
    int shearTag = 0;

    if (argc > 13) {
      for (int i = 13; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << endln;
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-constHinge") == 0 && ++i < argc) {
          if (Tcl_GetInt(interp, argv[i], &shearTag) != TCL_OK) {
            opserr << "WARNING invalid constHinge tag\n";
            opserr << "BeamWithHinges: " << tag << endln;
            return TCL_ERROR;
          }
          isShear = true;
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &numIters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            opserr << "BeamWithHinges: " << tag << endln;
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            opserr << "BeamWithHinges: " << tag << endln;
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    SectionForceDeformation *sectionI = builder->getTypedObject<SectionForceDeformation>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    SectionForceDeformation *sectionJ = builder->getTypedObject<SectionForceDeformation>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;

    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection2d elastic(8, E, A, I);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;

      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;
    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;
    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == nullptr) {
      opserr << "Unknown element type: " << argv[1] << endln;
      return TCL_ERROR;
    }

    if (isShear) {
      SectionForceDeformation *sectionL = builder->getTypedObject<SectionForceDeformation>(shearTag);
      if (sectionL == nullptr)
        return TCL_ERROR;

      sections[numSections++] = sectionL;
    }

    theElement = new ForceBeamColumn2d(tag, ndI, ndJ, numSections, sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       numIters, tol);

    delete theBeamIntegr;

    // Ensure we have created the element, out of memory if got here and no
    // element
    if (theElement == nullptr) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (theDomain->addElement(theElement) == false) {
      opserr << "WARNING TclElmtBuilder - addBeamWithHinges - could not add "
                "element to domain ";
      opserr << tag << endln;
      return TCL_ERROR;
    }
  }

  else if (NDM == 3 && NDF == 6) {
    if (argc < 16) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? Iz? Iy? G? J? transfTag? <-shear shearLength?> <-mass "
                "massDens?> <-iter maxIters tolerance>"
             << endln;
      return TCL_ERROR;
    }

    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, Iz, Iy, G, J;
    double massDens = 0.0;
    int numIters = 10;
    double tol = 1.0e-10;
    double shearLength = 1.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[12], &Iy) != TCL_OK) {
      opserr << "WARNING invalid Iy\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[13], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[14], &J) != TCL_OK) {
      opserr << "WARNING invalid J\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[15], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }


    if (argc > 16) {
      for (int i = 16; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << endln;
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-shear") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &shearLength) != TCL_OK) {
            opserr << "WARNING invalid shearLength\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &numIters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    SectionForceDeformation *sectionI = builder->getTypedObject<SectionForceDeformation>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    SectionForceDeformation *sectionJ = builder->getTypedObject<SectionForceDeformation>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;


    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection3d elastic(0, E, A, Iz, Iy, G, J);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;
    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;
    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == 0) {
      opserr << "Unknown element type: " << argv[1] << endln;
      return TCL_ERROR;
    }

    /*
    if (isShear) {
      SectionForceDeformation *sectionL = builder->getTypedObject<SectionForceDeformation>(shearTag);

      if (sectionL == 0) {
        opserr << "WARNING section L does not exist\n";
        opserr << "section: " << shearTag;
        opserr << "\nBeamWithHinges: " << tag << endln;
        return TCL_ERROR;
      }
      sections[numSections++] = sectionL;
    }
    */

    theElement = new ForceBeamColumn3d(tag, ndI, ndJ, numSections, sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       numIters, tol);

    delete theBeamIntegr;

    // Ensure we have created the element, out of memory if got here and no
    // element
    if (theElement == nullptr) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "BeamWithHinges: " << tag << endln;
      return TCL_ERROR;
    }

    if (theDomain->addElement(theElement) == false) {
      opserr << "WARNING TclElmtBuilder - addBeamWithHinges - could not add "
                "element to domain ";
      opserr << tag << endln;
      return TCL_ERROR;
    }
  }

  else {
    opserr << "ERROR -- model dimension: " << NDM
           << " and nodal degrees of freedom: " << NDF
           << " are incompatible for BeamWithHinges element" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}
