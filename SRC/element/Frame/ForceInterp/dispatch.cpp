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
//
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the
// TclBasicBuilder_addDispBeamColumn() command.
//
#include <utility>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

#include <assert.h>
#include <Domain.h>
#include <G3_Logging.h>

#include <runtime/BasicModelBuilder.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>
#include <ForceBeamColumnWarping2d.h>
#include <ElasticForceBeamColumnWarping2d.h>
#include <TimoshenkoBeamColumn2d.h>
#include <DispBeamColumn2d.h>
#include <DispBeamColumn3d.h>
#include <DispBeamColumnNL2d.h>
#include <DispBeamColumn2dThermal.h>
#include <DispBeamColumn3dThermal.h>  //L.Jiang [SIF]
#include <ForceBeamColumn2dThermal.h> //L.Jiang [SIF]

#include <CrdTransf.h>

#include <DispBeamColumn2dWithSensitivity.h>
#include <DispBeamColumn3dWithSensitivity.h>

#include <ElasticForceBeamColumn2d.h>
#include <ElasticForceBeamColumn3d.h>

#include <ForceBeamColumnCBDI2d.h>

#include <BeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <UserDefinedHingeIntegration.h>
#include <DistHingeIntegration.h>
#include <RegularizedHingeIntegration.h>

#include <TrapezoidalBeamIntegration.h>
#include <CompositeSimpsonBeamIntegration.h>
#include <FixedLocationBeamIntegration.h>
#include <LowOrderBeamIntegration.h>
#include <MidDistanceBeamIntegration.h>
//#include <GaussQBeamIntegration.h>

#include <ElasticSection2d.h>
#include <ElasticSection3d.h>


int
TclBasicBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,
                                   int inArgc, TCL_Char **const inArgv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain *domain = builder->getDomain();
  assert(domain != nullptr);

  bool deleteBeamIntegr = true;

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  int ok = 0;
  if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4))
    ok = 1;
  if (ndm == 3 && ndf == 6)
    ok = 1;

  if (ok == 0) {
    opserr << G3_ERROR_PROMPT << "ndm = " << ndm << " and ndf = " << ndf
           << " not compatible with Frame element" << endln;
    return TCL_ERROR;
  }

  // split possible lists present in argv
  char *List = Tcl_Merge(inArgc, inArgv);
  if (List == nullptr) {
    opserr << G3_ERROR_PROMPT << "problem merging list\n";
    return TCL_ERROR;
  }

  // remove braces from list
  for (int i = 0; List[i] != '\0'; i++) {
    if ((List[i] == '{') || (List[i] == '}'))
      List[i] = ' ';
  }

  int argc;
  TCL_Char ** argv;
  if (Tcl_SplitList(interp, List, &argc, &argv) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "problem splitting list\n";
    return TCL_ERROR;
  }
  Tcl_Free((char *)List);

  if (argc < 6) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
    opserr << "Want: element " << argv[1]
           << " eleTag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }

  int eleTag, iNode, jNode, transfTag;
  CrdTransf *theTransf2d      = nullptr;
  CrdTransf *theTransf3d      = nullptr;
  Element *theElement         = nullptr;
  BeamIntegration *beamIntegr = nullptr;

  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid " << " eleTag " << eleTag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid jNode\n";
    return TCL_ERROR;
  }

  //
  // fmk UNDOCUMENTED FEATURE -
  // all to take similar command to nonlinearBeamColumn & dispBeamColumn
  //
  // element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag 
  //     $transfTag <-mass $massDens> <-iter $maxIters $tol> 
  //     <-integration $intType>
  //
  if ((strcmp(argv[6], "Lobatto") != 0) && 
      (strcmp(argv[6], "Legendre") != 0) &&
      (strcmp(argv[6], "Radau") != 0) &&
      (strcmp(argv[6], "NewtonCotes") != 0) &&
      (strcmp(argv[6], "UserDefined") != 0) &&
      (strcmp(argv[6], "HingeMidpoint") != 0) &&
      (strcmp(argv[6], "HingeEndpoint") != 0) &&
      (strcmp(argv[6], "HingeRadau") != 0) &&
      (strcmp(argv[6], "HingeRadauTwo") != 0) &&
      (strcmp(argv[6], "UserHinge") != 0) &&
      (strcmp(argv[6], "DistHinge") != 0) &&
      (strcmp(argv[6], "RegularizedHinge") != 0) &&
      (strcmp(argv[6], "Trapezoidal") != 0) &&
      (strcmp(argv[6], "CompositeSimpson") != 0) &&
      (strcmp(argv[6], "FixedLocation") != 0) &&
      (strcmp(argv[6], "LowOrder") != 0) && 
      (strcmp(argv[6], "GaussQ") != 0) &&
      (strcmp(argv[6], "MidDistance") != 0)) {

    int nIP, secTag;

    if (Tcl_GetInt(interp, argv[5], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    if (nIP <= 0) {
      opserr << G3_ERROR_PROMPT << "invalid nIP, must be > 0\n";
      return TCL_ERROR;
    }

    int argi = 6;
    SectionForceDeformation **sections = nullptr;

    if (strcmp(argv[argi], "-sections") != 0) {

      if (Tcl_GetInt(interp, argv[argi], &secTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid secTag\n";
        return TCL_ERROR;
      } else
        argi++;

      // Check if we are being called from OpenSeesPy, in which case we need to parse things a little
      // differently
      const char *openseespy = Tcl_GetVar(interp, "opensees::pragma::openseespy", 0);
      if (openseespy == nullptr || strcmp(openseespy, "0")==0) { 

        // OpenSees Tcl behavior

        SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secTag);
        if (theSection == nullptr) {
          return TCL_ERROR;
        }
        sections = new SectionForceDeformation *[nIP];
        for (int i = 0; i < nIP; i++)
          sections[i] = theSection;

        // Geometric transformation
        if (Tcl_GetInt(interp, argv[argi], &transfTag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
          return TCL_ERROR;
        } else
          argi++;


      } else {

        // OpenSeesPy behavior with BeamIntegrations

        transfTag = nIP;
        BeamIntegrationRule* theRule = builder->getTypedObject<BeamIntegrationRule>(secTag);
        if (theRule == nullptr) {
          return TCL_ERROR;
        }

        deleteBeamIntegr = false;
        beamIntegr = theRule->getBeamIntegration();
        const ID& secTags = theRule->getSectionTags();
        nIP = secTags.Size();
        sections = new SectionForceDeformation *[nIP];
        for (int i=0; i < secTags.Size(); i++) {
          sections[i] = builder->getTypedObject<SectionForceDeformation>(secTags(i));
          if (sections[i] == nullptr) {
            opserr << "section " << secTags(i) << "not found\n";
            return TCL_ERROR;
          }
        }

      }

    } else {

      argi++;
      // get section tags
      sections = new SectionForceDeformation *[nIP];
      for (int i = 0; i < nIP; i++) {
        if (Tcl_GetInt(interp, argv[argi], &secTag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid secTag\n";
          return TCL_ERROR;
        } else
          argi++;

        SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secTag);
        if (theSection == nullptr) {
          opserr << G3_ERROR_PROMPT << "section not found\n";
          opserr << "Section: " << secTag;
          opserr << ", element: " << eleTag << endln;
          return TCL_ERROR;
        }

        sections[i] = theSection;
      }

      if (Tcl_GetInt(interp, argv[argi], &transfTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
        return TCL_ERROR;
      } else
        argi++;
    }

    int numIter = 10;
    double tol = 1.0e-12;
    double mass = 0.0;
    int cMass = 0;

    while (argi < argc) {
      if (strcmp(argv[argi], "-iter") == 0) {
        if (argc < argi + 3) {
          opserr << G3_ERROR_PROMPT << "not enough -iter args need -iter numIter? tol?\n";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &numIter) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid numIter\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 2], &tol) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid tol\n";
          return TCL_ERROR;
        }
        argi += 3;
      } else if (strcmp(argv[argi], "-mass") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough -mass args need -mass mass?\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid numIter\n";
          return TCL_ERROR;
        }
        argi += 2;
      } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                 (strcmp(argv[argi], "lMass") == 0)) {
        cMass = 0;
        argi++;
      } else if ((strcmp(argv[argi], "-cMass") == 0) ||
                 (strcmp(argv[argi], "cMass") == 0)) {
        cMass = 1;
        argi++;

      } else if (strcmp(argv[argi], "-integration") == 0) {

        argi++;
        if (strcmp(argv[argi], "Lobatto") == 0)
          beamIntegr = new LobattoBeamIntegration();
        else if (strcmp(argv[argi], "Legendre") == 0)
          beamIntegr = new LegendreBeamIntegration();
        else if (strcmp(argv[argi], "Radau") == 0)
          beamIntegr = new RadauBeamIntegration();
        else if (strcmp(argv[argi], "NewtonCotes") == 0)
          beamIntegr = new NewtonCotesBeamIntegration();
        else if (strcmp(argv[argi], "Trapezoidal") == 0)
          beamIntegr = new TrapezoidalBeamIntegration();
        else if (strcmp(argv[argi], "CompositeSimpson") == 0)
          beamIntegr = new CompositeSimpsonBeamIntegration();

        argi++;

        if (beamIntegr == nullptr) {
          opserr << G3_ERROR_PROMPT << "invalid integration type\n";
          return TCL_ERROR;
        }
      } else
        argi++;
    }

    //
    //
    //
    switch (ndm) {
      case 2:
        theTransf2d = builder->getTypedObject<CrdTransf>(transfTag);
        if (theTransf2d == nullptr) {
          opserr << G3_ERROR_PROMPT << "transformation not found with tag " << transfTag << endln;
          return TCL_ERROR;
        }
        break;

      case 3:
        theTransf3d = builder->getTypedObject<CrdTransf>(transfTag);
        if (theTransf3d == nullptr) {
          opserr << G3_ERROR_PROMPT << "transformation not found with tag " << transfTag << endln;
          return TCL_ERROR;
        }
    }

    if (beamIntegr == nullptr) {
      if (strstr(argv[1], "ispBeam") == 0) {
        beamIntegr = new LobattoBeamIntegration();
      } else {
        beamIntegr = new LegendreBeamIntegration();
      }
    }

    if (ndm == 2) {
      if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
        theElement = new ElasticForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "timoshenkoBeamColumn") == 0)
        theElement = new TimoshenkoBeamColumn2d  (eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "dispBeamColumn") == 0)
        theElement = new DispBeamColumn2d        (eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass, cMass);

      else if (strcmp(argv[1], "dispBeamColumnNL") == 0)
        theElement = new DispBeamColumnNL2d      (eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "forceBeamColumnCBDI") == 0)
        theElement = new ForceBeamColumnCBDI2d   (eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "forceBeamColumnCSBDI") == 0)
        theElement = new ForceBeamColumnCBDI2d   (eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass, true);

      else if (strcmp(argv[1], "forceBeamColumnWarping") == 0)
        theElement = new ForceBeamColumnWarping2d(eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "elasticForceBeamColumnWarping") == 0)
        theElement = new ElasticForceBeamColumnWarping2d(eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);
      else if (strcmp(argv[1], "dispBeamColumnThermal") == 0)
        theElement =
            new DispBeamColumn2dThermal(eleTag, iNode, jNode, nIP, sections,
                                        *beamIntegr, *theTransf2d, mass);
      else if (strcmp(argv[1], "forceBeamColumnThermal") == 0)
        theElement = new ForceBeamColumn2dThermal(eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "dispBeamColumnWithSensitivity") == 0)
        theElement = new DispBeamColumn2dWithSensitivity(
            eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf2d,
            mass);
      else
        theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
                                           *beamIntegr, *theTransf2d, mass,
                                           numIter, tol);
    } else {
      if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
        theElement =
            new ElasticForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
                                         *beamIntegr, *theTransf3d, mass);
      else if (strcmp(argv[1], "dispBeamColumn") == 0)
        theElement =
            new DispBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
                                 *beamIntegr, *theTransf3d, mass, cMass);
      else if (strcmp(argv[1], "dispBeamColumnThermal") == 0)
        theElement = new DispBeamColumn3dThermal(
            eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf3d, mass);
      else if (strcmp(argv[1], "dispBeamColumnWithSensitivity") == 0)
        theElement = new DispBeamColumn3dWithSensitivity(
            eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf3d, mass);
      else
        theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
                                           *beamIntegr, *theTransf3d, mass,
                                           numIter, tol);
    }

    delete[] sections;
    if (deleteBeamIntegr)
      delete beamIntegr;

    if (domain->addElement(theElement) == false) {
      opserr << G3_ERROR_PROMPT << "could not add element to the domain\n";
      delete theElement;
      return TCL_ERROR;
    }

    Tcl_Free((char *)argv);

    return TCL_OK;
  }

  //
  // otherwise use correct format of command as found in current documentation
  //
  //  element forceBeamColumn $eleTag $iNode $jNode 
  //          $transfTag "IntegrationType arg1 arg2 ..." 
  //          <-mass $massDens> <-iter $maxIters $tol>
  //
  if (Tcl_GetInt(interp, argv[5], &transfTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
    return TCL_ERROR;
  }

  if (ndm == 2) {
    theTransf2d = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf2d == nullptr) {
      opserr << G3_ERROR_PROMPT << "transformation not found\n";
      return TCL_ERROR;
    }
  } else if (ndm == 3) {
    theTransf3d = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf3d == nullptr) {
      opserr << G3_ERROR_PROMPT << "transformation not found\n";
      return TCL_ERROR;
    }
  }

//BeamIntegration *beamIntegr = nullptr;
  SectionForceDeformation **sections;
  int numSections;

  if (strcmp(argv[6], "Lobatto") == 0 || strcmp(argv[6], "Legendre") == 0 ||
      strcmp(argv[6], "Radau") == 0 || strcmp(argv[6], "NewtonCotes") == 0 ||
      strcmp(argv[6], "Trapezoidal") == 0 ||
      strcmp(argv[6], "CompositeSimpson") == 0) {

    int secTag;

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secTag);
    if (theSection == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found with tag " << secTag << endln;
      return TCL_ERROR;
    }

    sections = new SectionForceDeformation *[numSections];
    for (int i = 0; i < numSections; i++)
      sections[i] = theSection;

    if (strcmp(argv[6], "Lobatto") == 0)
      beamIntegr = new LobattoBeamIntegration();
    else if (strcmp(argv[6], "Legendre") == 0)
      beamIntegr = new LegendreBeamIntegration();
    else if (strcmp(argv[6], "Radau") == 0)
      beamIntegr = new RadauBeamIntegration();
    else if (strcmp(argv[6], "NewtonCotes") == 0)
      beamIntegr = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[6], "Trapezoidal") == 0)
      beamIntegr = new TrapezoidalBeamIntegration();
    else if (strcmp(argv[6], "CompositeSimpson") == 0)
      beamIntegr = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[6] << endln;
      return TCL_ERROR;
    }
  }
#if 0
  else if (strcmp(argv[6],"GaussQ") == 0) {

    int type, secTag;

    if (argc < 10) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1] << " eleTag? iNode? jNode? transfTag? GaussQ type? secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &type) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid type\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &secTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secTag);
    if (theSection == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTag;
      return TCL_ERROR;
    }

    sections = new SectionForceDeformation *[numSections];
    for (int i = 0; i < numSections; i++)
      sections[i] = theSection;

    beamIntegr = new GaussQBeamIntegration(type);
  }
#endif
  else if (strcmp(argv[6], "UserDefined") == 0) {

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? UserDefined nIP? secTag1? "
                "... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSections\n";
      return TCL_ERROR;
    }

    ID secs(numSections);
    Vector pts(numSections);
    Vector wts(numSections);

    int i, j;
    for (i = 0, j = 8; i < numSections; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + numSections], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * numSections], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i) = pt;
      wts(i) = wt;
    }

    sections = new SectionForceDeformation *[numSections];
    for (i = 0; i < numSections; i++) {
      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secs(i));
      if (theSection == nullptr) {
        opserr << G3_ERROR_PROMPT << "section not found\n";
        opserr << "Section: " << secs(i);
        return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    beamIntegr = new UserDefinedBeamIntegration(numSections, pts, wts);
  }

  else if (strcmp(argv[6], "HingeMidpoint") == 0 ||
           strcmp(argv[6], "HingeRadau") == 0 ||
           strcmp(argv[6], "HingeRadauTwo") == 0 ||
           strcmp(argv[6], "HingeEndpoint") == 0) {

    if (argc < 12) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? type secTagI? lpI? secTagJ? "
                "lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;

    if (Tcl_GetInt(interp, argv[7], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[11], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionI = builder->getTypedObject<SectionForceDeformation>(secTagI);
    if (sectionI == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagI;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = builder->getTypedObject<SectionForceDeformation>(secTagJ);
    if (sectionJ == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagJ;
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionE = builder->getTypedObject<SectionForceDeformation>(secTagE);
    if (sectionJ == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagE;
      return TCL_ERROR;
    }

    sections = new SectionForceDeformation *[6];

    if (strcmp(argv[6], "HingeMidpoint") == 0) {
      beamIntegr = new HingeMidpointBeamIntegration(lpI, lpJ);
      numSections = 4;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionJ;
    } else if (strcmp(argv[6], "HingeRadau") == 0) {
      beamIntegr = new HingeRadauBeamIntegration(lpI, lpJ);
      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionE;
      sections[4] = sectionE;
      sections[5] = sectionJ;
    } else if (strcmp(argv[6], "HingeRadauTwo") == 0) {
      beamIntegr = new HingeRadauTwoBeamIntegration(lpI, lpJ);
      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = sectionE;
      sections[3] = sectionE;
      sections[4] = sectionJ;
      sections[5] = sectionJ;
    } else {
      beamIntegr = new HingeEndpointBeamIntegration(lpI, lpJ);
      numSections = 4;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionJ;
    }
  }

  else if (strcmp(argv[6], "UserHinge") == 0) {

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? UserHinge secTagE? npL? "
                "secTagL1? ... ptL1? ... wtL1? ... npR? secTagR1? ... ptR1? "
                "... wtR1? ...\n";
      return TCL_ERROR;
    }

    int secTagE;

    if (Tcl_GetInt(interp, argv[7], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    int argStart = 8;

    int npL, npR;

    if (Tcl_GetInt(interp, argv[argStart], &npL) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid npL\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart + 3 * npL + 1], &npR) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid npR\n";
      return TCL_ERROR;
    }

    numSections = npL + npR;

    ID secs(numSections);
    Vector ptsL(npL);
    Vector wtsL(npL);
    Vector ptsR(npR);
    Vector wtsR(npR);

    int i, j;
    for (i = 0, j = argStart + 1; i < npL; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + npL], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * npL], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      secs(i) = sec;
      ptsL(i) = pt;
      wtsL(i) = wt;
    }
    for (i = 0, j = 1 + (argStart + 1) + 3 * npL; i < npR; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + npR], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * npR], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      secs(i + npL) = sec;
      ptsR(i) = pt;
      wtsR(i) = wt;
    }

    sections = new SectionForceDeformation *[numSections + 2];

    for (i = 0; i < numSections; i++) {
      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secs(i));
      if (theSection == nullptr) {
        opserr << G3_ERROR_PROMPT << "section not found\n";
        opserr << "Section: " << secs(i);
        return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    SectionForceDeformation *sectionE = builder->getTypedObject<SectionForceDeformation>(secTagE);
    if (sectionE == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagE;
      return TCL_ERROR;
    }

    sections[numSections] = sectionE;
    sections[numSections + 1] = sectionE;

    beamIntegr =
        new UserDefinedHingeIntegration(npL, ptsL, wtsL, npR, ptsR, wtsR);

    numSections += 2;
  }

  else if (strcmp(argv[6], "DistHinge") == 0) {

    if (argc < 14) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? type distType nIP? secTagI? "
                "lpI? secTagJ? lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;
    int nIP;

    BeamIntegration *otherBeamInt = 0;
    if (strcmp(argv[7], "Lobatto") == 0)
      otherBeamInt = new LobattoBeamIntegration();
    else if (strcmp(argv[7], "Legendre") == 0)
      otherBeamInt = new LegendreBeamIntegration();
    else if (strcmp(argv[7], "Radau") == 0)
      otherBeamInt = new RadauBeamIntegration();
    else if (strcmp(argv[7], "NewtonCotes") == 0)
      otherBeamInt = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[7], "Trapezoidal") == 0)
      otherBeamInt = new TrapezoidalBeamIntegration();
    else if (strcmp(argv[7], "CompositeSimpson") == 0)
      otherBeamInt = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[7] << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[11], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[13], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionI = builder->getTypedObject<SectionForceDeformation>(secTagI);
    if (sectionI == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagI;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = builder->getTypedObject<SectionForceDeformation>(secTagJ);
    if (sectionJ == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagJ;
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionE = builder->getTypedObject<SectionForceDeformation>(secTagE);
    if (sectionJ == nullptr) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagE;
      return TCL_ERROR;
    }

    numSections = 2 * nIP;
    sections = new SectionForceDeformation *[numSections + 2];
    for (int i = 0; i < nIP; i++) {
      sections[i] = sectionI;
      sections[i + nIP] = sectionJ;
    }

    sections[numSections] = sectionE;
    sections[numSections + 1] = sectionE;

    beamIntegr = new DistHingeIntegration(lpI, lpJ, *otherBeamInt);

    numSections += 2;

    if (otherBeamInt != 0)
      delete otherBeamInt;
  }

  else if (strcmp(argv[6], "RegularizedHinge") == 0) {

    if (argc < 16) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? type distType nIP? secTagI? "
                "lpI? zetaI? secTagJ? lpJ? zetaJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;
    double zetaI, zetaJ;
    int nIP;

    BeamIntegration *otherBeamInt = nullptr;
    if (strcmp(argv[7], "Lobatto") == 0)
      otherBeamInt = new LobattoBeamIntegration();
    else if (strcmp(argv[7], "Legendre") == 0)
      otherBeamInt = new LegendreBeamIntegration();
    else if (strcmp(argv[7], "Radau") == 0)
      otherBeamInt = new RadauBeamIntegration();
    else if (strcmp(argv[7], "NewtonCotes") == 0)
      otherBeamInt = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[7], "Trapezoidal") == 0)
      otherBeamInt = new TrapezoidalBeamIntegration();
    else if (strcmp(argv[7], "CompositeSimpson") == 0)
      otherBeamInt = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[7] << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &zetaI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zetaI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[12], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[14], &zetaJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zetaI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[15], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionI = builder->getTypedObject<SectionForceDeformation>(secTagI);
    if (sectionI == 0) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagI;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = builder->getTypedObject<SectionForceDeformation>(secTagJ);
    if (sectionJ == 0) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagJ;
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionE = builder->getTypedObject<SectionForceDeformation>(secTagE);
    if (sectionJ == 0) {
      opserr << G3_ERROR_PROMPT << "section not found\n";
      opserr << "Section: " << secTagE;
      return TCL_ERROR;
    }

    numSections = nIP;
    sections = new SectionForceDeformation *[numSections + 2];
    for (int i = 1; i < nIP - 1; i++) {
      sections[i] = sectionE;
    }

    sections[0] = sectionI;
    sections[numSections] = sectionI;
    sections[numSections - 1] = sectionJ;
    sections[numSections + 1] = sectionJ;

    beamIntegr =
        new RegularizedHingeIntegration(*otherBeamInt, lpI, lpJ, zetaI, zetaJ);

    numSections += 2;

    if (otherBeamInt != 0)
      delete otherBeamInt;
  }

  else if (strcmp(argv[6], "FixedLocation") == 0) {

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? FixedLocation nIP? secTag1? "
                "... pt1? ... \n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSections\n";
      return TCL_ERROR;
    }

    ID secs(numSections);
    Vector pts(numSections);
    Vector wts(numSections);

    int i, j;
    for (i = 0, j = 8; i < numSections; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + numSections], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i) = pt;
    }

    sections = new SectionForceDeformation *[numSections];
    for (i = 0; i < numSections; i++) {
      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secs(i));
      if (theSection == nullptr) {
        opserr << G3_ERROR_PROMPT << "section not found\n";
        opserr << "Section: " << secs(i);
        return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    beamIntegr = new FixedLocationBeamIntegration(numSections, pts);
  }

  else if (strcmp(argv[6], "LowOrder") == 0) {

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? LowOrder nIP? secTag1? ... "
                "pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid numSections\n";
      return TCL_ERROR;
    }

    ID secs(numSections);
    Vector pts(numSections);
    Vector wts(numSections);

    int i, j;
    int nc = 0;
    for (i = 0, j = 8; i < numSections; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      secs(i) = sec;

      if (Tcl_GetDouble(interp, argv[j + numSections], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      pts(i) = pt;

      if (j + 2 * numSections < argc) {
        if (Tcl_GetDouble(interp, argv[j + 2 * numSections], &wt) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid wt\n";
          return TCL_ERROR;
        } else {
          wts(i) = wt;
          nc++;
        }
      }
    }

    sections = new SectionForceDeformation *[numSections];
    for (i = 0; i < numSections; i++) {
      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secs(i));
      if (theSection == nullptr) {
        opserr << G3_ERROR_PROMPT << "section not found\n";
        opserr << "Section: " << secs(i);
        return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    beamIntegr = new LowOrderBeamIntegration(numSections, pts, nc, wts);
  }

  else if (strcmp(argv[6], "MidDistance") == 0) {

    if (argc < 9) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: element " << argv[1]
             << " eleTag? iNode? jNode? transfTag? MidDistance nIP? secTag1? "
                "... pt1? ... \n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &numSections) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    ID secs(numSections);
    Vector pts(numSections);
    Vector wts(numSections);

    int i, j;
    for (i = 0, j = 8; i < numSections; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + numSections], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i) = pt;
    }

    sections = new SectionForceDeformation *[numSections];
    for (i = 0; i < numSections; i++) {
      SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secs(i));
      if (theSection == nullptr) {
        opserr << G3_ERROR_PROMPT << "section not found\n";
        opserr << "Section: " << secs(i);
        return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    // Sort locations before calling constructor
    for (int i = 0; i < numSections; i++) {
      int key = i;
      for (int j = i + 1; j < numSections; j++) {
        if (pts(j) < pts(key)) {
          key = j;
        }
      }
      if (key != i) {
        // Swap locs
        double temp;
        temp = pts(i);
        pts(i) = pts(key);
        pts(key) = temp;
        // Swap sections
        SectionForceDeformation *tempSection;
        tempSection = sections[i];
        sections[i] = sections[key];
        sections[key] = tempSection;
      }
    }

    beamIntegr = new MidDistanceBeamIntegration(numSections, pts);
  }

  else {
    opserr << "Unknown integration type: " << argv[6] << endln;
    opserr << argv[1] << " element: " << eleTag << endln;
    return TCL_ERROR;
  }

  int argi = 6;
  double mass = 0.0;
  int cMass = 0;
  int numIter = 10;
  double tol = 1.0e-12;

  while (argi < argc) {
    if (strcmp(argv[argi], "-iter") == 0) {
      if (argc < argi + 3) {
        opserr << G3_ERROR_PROMPT << "not enough -iter args need -iter numIter? tol?\n";
        opserr << argv[1] << " element: " << eleTag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[argi + 1], &numIter) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid numIter\n";
        opserr << argv[1] << " element: " << eleTag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[argi + 2], &tol) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid tol\n";
        opserr << argv[1] << " element: " << eleTag << endln;
        return TCL_ERROR;
      }
      argi += 2;
    } else if (strcmp(argv[argi], "-mass") == 0) {
      if (argc < argi + 2) {
        opserr << G3_ERROR_PROMPT << "not enough -mass args need -mass mass?\n";
        opserr << argv[1] << " element: " << eleTag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid numIter\n";
        opserr << argv[1] << " element: " << eleTag << endln;
        return TCL_ERROR;
      }
      argi += 1;
    } else if ((strcmp(argv[argi], "-lMass") == 0 ||
                strcmp(argv[argi], "lMass") == 0)) {
      cMass = 0;
    } else if ((strcmp(argv[argi], "-cMass") == 0 ||
                strcmp(argv[argi], "cMass") == 0)) {
      cMass = 1;
    }
    argi += 1;
  }

  if (ndm == 2) {
    if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn2d(
          eleTag, iNode, jNode, numSections, sections, *beamIntegr,
          *theTransf2d, mass);
    else if (strcmp(argv[1], "timoshenkoBeamColumn") == 0)
      theElement =
          new TimoshenkoBeamColumn2d(eleTag, iNode, jNode, numSections,
                                     sections, *beamIntegr, *theTransf2d, mass);
    else if (strcmp(argv[1], "dispBeamColumn") == 0)
      theElement =
          new DispBeamColumn2d(eleTag, iNode, jNode, numSections, sections,
                               *beamIntegr, *theTransf2d, mass, cMass);
    else if (strcmp(argv[1], "dispBeamColumnNL") == 0)
      theElement =
          new DispBeamColumnNL2d(eleTag, iNode, jNode, numSections, sections,
                                 *beamIntegr, *theTransf2d, mass);
    else if (strcmp(argv[1], "forceBeamColumnCBDI") == 0)
      theElement =
          new ForceBeamColumnCBDI2d(eleTag, iNode, jNode, numSections, sections,
                                    *beamIntegr, *theTransf2d, mass);
    else if (strcmp(argv[1], "forceBeamColumnCSBDI") == 0)
      theElement =
          new ForceBeamColumnCBDI2d(eleTag, iNode, jNode, numSections, sections,
                                    *beamIntegr, *theTransf2d, mass, true);
    else if (strcmp(argv[1], "forceBeamColumnWarping") == 0)
      theElement =
          new ForceBeamColumnWarping2d(eleTag, iNode, jNode, numSections,
                                       sections, *beamIntegr, *theTransf2d);
    else if (strcmp(argv[1], "elasticForceBeamColumnWarping") == 0)
      theElement = new ElasticForceBeamColumnWarping2d(
          eleTag, iNode, jNode, numSections, sections, *beamIntegr,
          *theTransf2d);
    else
      theElement =
          new ForceBeamColumn2d(eleTag, iNode, jNode, numSections, sections,
                                *beamIntegr, *theTransf2d, mass, numIter, tol);
  } else {
    if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn3d(
          eleTag, iNode, jNode, numSections, sections, *beamIntegr,
          *theTransf3d, mass);

    else if (strcmp(argv[1], "dispBeamColumn") == 0)
      theElement =
          new DispBeamColumn3d(eleTag, iNode, jNode, numSections, sections,
                               *beamIntegr, *theTransf3d, mass, cMass);
    else
      theElement =
          new ForceBeamColumn3d(eleTag, iNode, jNode, numSections, sections,
                                *beamIntegr, *theTransf3d, mass, numIter, tol);
  }

  if (beamIntegr != nullptr)
    delete beamIntegr;

  if (sections != nullptr)
    delete[] sections;

  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "ran out of memory creating element\n";
    return TCL_ERROR;
  }

  if (domain->addElement(theElement) == false) {
    opserr << G3_ERROR_PROMPT << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  Tcl_Free((char *)argv);

  return TCL_OK;
}
