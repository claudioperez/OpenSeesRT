/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the function to parse the TCL input
//              for the ElasticBeamColumn element.
//
// Written: fmk
// Created: 07/99
//
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <CrdTransf.h>
#include "ElasticBeam2d.h"
#include "ElasticBeam3d.h"
#include "PrismFrame3d.h"
#include <SectionForceDeformation.h>

#include <runtime/BasicModelBuilder.h>


int
TclBasicBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{

  // ensure the destructor has not been called
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  struct BeamData {
    double E, G, A, Iz, Iy, Ixy, J, Cw;
  } beamData;
  int beamId, iNode, jNode, transTag;

  Element *theBeam = nullptr;

  int eleArgStart = 1;
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  //
  // Preliminary checks
  //

  // Check problem dimension
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING elasticBeamColumn command only works when ndm is 2 or "
              "3, ndm: ";
    opserr << ndm << "\n";
    return TCL_ERROR;
  }
  // Check the number of arguments
  if ((argc - eleArgStart) < 4) {
    opserr << "WARNING bad command - want: element " << argv[1] << " beamId iNode jNode ...\n";
    return TCL_ERROR;
  }

  //
  // Parse common parameters
  //

  // tag, iNode, jNode, E
  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &beamId) != TCL_OK) {
    opserr << "WARNING invalid beamId: " << argv[1 + eleArgStart];
    opserr << " - elasticBeamColumn beamId iNode jNode A E I\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode - elasticBeamColumn " << beamId
           << " iNode jNode A E I\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode - elasticBeamColumn " << beamId
           << " iNode jNode A E I\n";
    return TCL_ERROR;
  }

  //
  // Parse unique to 2D or 3D
  //

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - ElasticBeamColumn \n";
      return TCL_ERROR;
    }

    // check the number of arguments
    if ((argc - eleArgStart) < 8) {
      opserr << "WARNING bad command - want: ElasticBeamColumn beamId iNode "
                "jNode A E I <alpha> <d> transTag <-mass m> <-cMass>\n";
      return TCL_ERROR;
    }

    // get the section properties
    double A, E, I;
    if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &A) != TCL_OK) {
      opserr << "WARNING invalid A - elasticBeamColumn " << beamId
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &E) != TCL_OK) {
      opserr << "WARNING invalid E - elasticBeam " << beamId
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &I) != TCL_OK) {
      opserr << "WARNING invalid I - elasticBeamColumn " << beamId
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    double alpha = 0.0;
    double d = 0.0;
    int argi;
    double mass = 0.0;
    int cMass = 0;
    if (((argc - eleArgStart) == 10) &&
        (strcmp(argv[eleArgStart + 8], "-mass") != 0)) {

      if (Tcl_GetDouble(interp, argv[7 + eleArgStart], &alpha) != TCL_OK) {
        opserr << "WARNING invalid alpha - ElasticBeamColumn " << beamId
               << " iNode jNode A E I alpha depth \n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &d) != TCL_OK) {
        opserr << "WARNING invalid d - elasticBeamColumn " << beamId
               << " iNode jNode A E I alpha d \n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[9 + eleArgStart], &transTag) != TCL_OK) {
        opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId
               << " iNode jNode A E I alpha d transTag\n";
        return TCL_ERROR;
      }
      argi = 10;
    }

    else {
      if (Tcl_GetInt(interp, argv[7 + eleArgStart], &transTag) != TCL_OK) {
        opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId
               << " iNode jNode A E I transTag\n";
        return TCL_ERROR;
      }
      argi = 8;
    }

    CrdTransf *theTrans = builder->getTypedObject<CrdTransf>(transTag);
    if (theTrans == nullptr)
      return TCL_ERROR;

    while (argi < argc) {
      if (strcmp(argv[argi], "-mass") == 0) {
        if (argc < argi + 2) {
          opserr << "WARNING not enough -mass args need -mass mass??\n";
          opserr << argv[1] << " element: " << beamId << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid mass\n";
          opserr << argv[1] << " element: " << beamId << endln;
          return TCL_ERROR;
        }
        argi += 2;
      } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                 (strcmp(argv[argi], "lMass") == 0)) {
        cMass = 0; // lumped mass matrix (default)
        argi++;
      } else if ((strcmp(argv[argi], "-cMass") == 0) ||
                 (strcmp(argv[argi], "cMass") == 0)) {
        cMass = 1; // consistent mass matrix
        argi++;
      } else
        argi++;
    }

    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam2d(beamId, A, E, I, iNode, jNode, *theTrans, alpha,
                                d, mass, cMass);
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndof: " << ndf;
      opserr << ", for 3d problem  need 6 - elasticBeamColumn \n";
      return TCL_ERROR;
    }

    // check the number of arguments
    if (((argc - eleArgStart) < 11) && ((argc - eleArgStart) != 6)) {
      opserr
          << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode";
      opserr << " A E G Jx Iy Iz transTag <-mass m> <-cMass>" << endln;
      return TCL_ERROR;
    }

    // get the section properties
    double A, E, G, Jx, Iy, Iz;
    if ((argc - eleArgStart) == 6) {
      int section;
      if (Tcl_GetInt(interp, argv[4 + eleArgStart], &section) != TCL_OK) {
        opserr << "WARNING invalid secTag - elasticBeamColumn iNode jNode "
                  "sectionTag? transfTag?\n";
        return TCL_ERROR;
      }

      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(section);
      if (theSection == nullptr)
        return TCL_ERROR;

      if (Tcl_GetInt(interp, argv[5 + eleArgStart], &transTag) != TCL_OK) {
        opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
        opserr << " iNode jNode sectionTag? transfTag?\n";
        return TCL_ERROR;
      }

      double mass = 0.0;
      int cMass = 0;
      int argi = 6 + eleArgStart;

      while (argi < argc) {
        if (strcmp(argv[argi], "-mass") == 0) {
          if (argc < argi + 2) {
            opserr << "WARNING not enough -mass args need -mass mass??\n";
            opserr << argv[1] << " element: " << beamId << endln;
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
            opserr << "WARNING invalid mass\n";
            opserr << argv[1] << " element: " << beamId << endln;
            return TCL_ERROR;
          }
          argi += 2;
        } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                   (strcmp(argv[argi], "lMass") == 0)) {
          cMass = 0; // lumped mass matrix (default)
          argi++;
        } else if ((strcmp(argv[argi], "-cMass") == 0) ||
                   (strcmp(argv[argi], "cMass") == 0)) {
          cMass = 1; // consistent mass matrix
          argi++;
        } else
          argi++;
      }


      CrdTransf *theTrans = builder->getTypedObject<CrdTransf>(transTag);
      if (theTrans == nullptr)
        return TCL_ERROR;

      // now create the beam and add it to the Domain
//    theBeam = new ElasticBeam3d(beamId, iNode, jNode, theSection, *theTrans, mass, cMass);
      theBeam = new PrismFrame3d(beamId, iNode, jNode, theSection, *theTrans, mass, cMass);

    } else {

      if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &A) != TCL_OK) {
        opserr << "WARNING invalid A - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &E) != TCL_OK) {
        opserr << "WARNING invalid E - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &G) != TCL_OK) {
        opserr << "WARNING invalid G - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[7 + eleArgStart], &Jx) != TCL_OK) {
        opserr << "WARNING invalid Jx - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &Iy) != TCL_OK) {
        opserr << "WARNING invalid Iy - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[9 + eleArgStart], &Iz) != TCL_OK) {
        opserr << "WARNING invalid Iz - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[10 + eleArgStart], &transTag) != TCL_OK) {
        opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
        opserr << " iNode jNode A E G Jx Iy Iz\n";
        return TCL_ERROR;
      }

      double mass = 0.0;
      int cMass = 0;
      int argi = 11 + eleArgStart;

      while (argi < argc) {
        if (strcmp(argv[argi], "-mass") == 0) {
          if (argc < argi + 2) {
            opserr << "WARNING not enough -mass args need -mass mass??\n";
            opserr << argv[1] << " element: " << beamId << endln;
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
            opserr << "WARNING invalid mass\n";
            opserr << argv[1] << " element: " << beamId << endln;
            return TCL_ERROR;
          }
          argi += 2;
        } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                   (strcmp(argv[argi], "lMass") == 0)) {
          cMass = 0; // lumped mass matrix (default)
          argi++;
        } else if ((strcmp(argv[argi], "-cMass") == 0) ||
                   (strcmp(argv[argi], "cMass") == 0)) {
          cMass = 1; // consistent mass matrix
          argi++;
        } else
          argi++;
      }

      CrdTransf *theTrans = builder->getTypedObject<CrdTransf>(transTag);
      if (theTrans == nullptr)
        return TCL_ERROR;


      // now create the beam and add it to the Domain
      theBeam = new ElasticBeam3d(beamId, A, E, G, Jx, Iy, Iz, iNode, jNode,
                                  *theTrans, mass, cMass);
//    theBeam = new PrismFrame3d (beamId, A, E, G, Jx, Iy, Iz, iNode, jNode,
//                                *theTrans, mass, cMass);
    }
  }


  // now add the beam to the domain
  if (builder->getDomain()->addElement(theBeam) == false) {
    opserr
        << "WARNING TclBasicBuilder - addBeam - could not add beam to domain ";
    opserr << *theBeam;
    delete theBeam;
    return TCL_ERROR;
  }

  return TCL_OK;
}
