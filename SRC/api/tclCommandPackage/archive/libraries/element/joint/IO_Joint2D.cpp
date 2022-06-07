

#include <g3_api.h>


#include <SRC/element/joint/Joint2D.h>
void *OPS_Joint2D()
{
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata != 8 && numdata != 10 && numdata != 12 && numdata != 18) {
    opserr << "WARNING incorrect number of arguments\n";
    opserr << "Want:\n";
    opserr
        << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? "
              "LrgDsp? -damage DmgTag?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? "
              "MatK? MatL? MatC? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? "
              "MatK? MatL? MatC? LrgDsp? -damage DmgI DmgJ DmgK DmgL DmgC\n";
    return 0;
  }

  // Joint2DId, iNode, jNode, kNode, lNode, CenterNodeTag
  int idata[6];
  int num = 6;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data\n";
    return 0;
  }
  int Joint2DId = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];
  int kNode = idata[3];
  int lNode = idata[4];
  int CenterNodeTag = idata[5];

  // check domain for existence of internal node tag
  Node *CenterNode = theDomain->getNode(CenterNodeTag);
  if (CenterNode != 0) {
    opserr
        << "WARNING node tag specified for the center node already exists.\n";
    opserr << "Use a new node tag.\n";
    opserr << "Joint2D element: " << Joint2DId << endln;
    return 0;
  }

  UniaxialMaterial *MatI = NULL;
  UniaxialMaterial *MatJ = NULL;
  UniaxialMaterial *MatK = NULL;
  UniaxialMaterial *MatL = NULL;
  UniaxialMaterial *PanelMaterial = NULL;
  Joint2D *theJoint2D;
  int LargeDisp;

  // Decide to use which constructor, based on the number of arguments
  if (numdata == 8 || numdata == 12) {

    // Using Joint2D constructor without damage

    if (numdata == 8) {
      int PanelMatId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelMatId) < 0) {
        opserr << "WARNING invalid matID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
        // use 0 as default
        LargeDisp = 0;
      }

      PanelMaterial = OPS_getUniaxialMaterial(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }
    }

    else // if ( (argc-argStart) == 12  )
    {
      int MatIid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatIid) < 0) {
        opserr << "WARNING invalid material ID for spring I\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (MatIid != 0) {
        MatI = OPS_getUniaxialMaterial(MatIid);

        if (MatI == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatIid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatI = NULL;

      int MatJid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatJid) < 0) {
        opserr << "WARNING invalid material ID for spring J\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (MatJid != 0) {
        MatJ = OPS_getUniaxialMaterial(MatJid);

        if (MatJ == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatJid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatJ = NULL;

      int MatKid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatKid) < 0) {
        opserr << "WARNING invalid material ID for spring K\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }
      if (MatKid != 0) {
        MatK = OPS_getUniaxialMaterial(MatKid);

        if (MatK == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatKid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatK = NULL;

      int MatLid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatLid) < 0) {
        opserr << "WARNING invalid material ID for spring L\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }
      if (MatLid != 0) {
        MatL = OPS_getUniaxialMaterial(MatLid);

        if (MatL == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatL = NULL;

      int PanelMatId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelMatId) < 0) {
        opserr << "WARNING invalid matID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }
      PanelMaterial = OPS_getUniaxialMaterial(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }

      num = 1;
      if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
        // use 0 as default
        LargeDisp = 0;
      }
    }

    UniaxialMaterial *springModels[5] = {MatI, MatJ, MatK, MatL, PanelMaterial};
    theJoint2D = new Joint2D(Joint2DId, iNode, jNode, kNode, lNode,
                             CenterNodeTag, springModels, theDomain, LargeDisp);

    return theJoint2D;

  }

  else if (numdata == 10 || numdata == 18) {
    // Using Joint2D constructor with damage
    DamageModel *DmgI = NULL;
    DamageModel *DmgJ = NULL;
    DamageModel *DmgK = NULL;
    DamageModel *DmgL = NULL;
    DamageModel *PanelDamage = NULL;

    if (numdata == 10) {
      int PanelMatId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelMatId) < 0) {
        opserr << "WARNING invalid matID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      num = 1;
      if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
        // use 0 as default
        LargeDisp = 0;
      }

      PanelMaterial = OPS_getUniaxialMaterial(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }

      const char *damageFlag = OPS_GetString();
      if (strcmp(damageFlag, "-damage") != 0 &&
          strcmp(damageFlag, "-Damage") != 0) {
        opserr << "WARNING incorrect command line\n";
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }

      int PanelDamageId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelDamageId) < 0) {
        opserr << "WARNING invalid damageID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      DamageModel *PanelDamage;
      PanelDamage = OPS_getDamageModel(PanelDamageId);

      if (PanelDamage == 0) {
        opserr << "WARNING damage model not found\n";
        opserr << "Damage Model: " << PanelDamageId;
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }
    }

    else // if ( (argc-argStart) == 18  )
    {
      int MatIid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatIid) < 0) {
        opserr << "WARNING invalid material ID for spring I\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (MatIid != 0) {
        MatI = OPS_getUniaxialMaterial(MatIid);

        if (MatI == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatIid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatI = NULL;

      int MatJid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatJid) < 0) {
        opserr << "WARNING invalid material ID for spring J\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (MatJid != 0) {
        MatJ = OPS_getUniaxialMaterial(MatJid);

        if (MatJ == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatJid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatJ = NULL;

      int MatKid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatKid) < 0) {
        opserr << "WARNING invalid material ID for spring K\n";
        opserr << "Joint2D element: " << Joint2DId << endln;

        return 0;
      }
      if (MatKid != 0) {
        MatK = OPS_getUniaxialMaterial(MatKid);

        if (MatK == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatKid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatK = NULL;

      int MatLid;
      num = 1;
      if (OPS_GetIntInput(&num, &MatLid) < 0) {
        opserr << "WARNING invalid material ID for spring L\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }
      if (MatLid != 0) {
        MatL = OPS_getUniaxialMaterial(MatLid);

        if (MatL == NULL) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        MatL = NULL;

      int PanelMatId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelMatId) < 0) {
        opserr << "WARNING invalid matID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }
      PanelMaterial = OPS_getUniaxialMaterial(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }

      num = 1;
      if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
        // use 0 as default
        LargeDisp = 0;
      }

      const char *damageFlag = OPS_GetString();
      if (strcmp(damageFlag, "-damage") != 0 &&
          strcmp(damageFlag, "-Damage") != 0) {
        opserr << "WARNING incorrect command line\n";
        opserr << "\nJoint2D element: " << Joint2DId << endln;
        return 0;
      }

      int DmgIid;
      num = 1;
      if (OPS_GetIntInput(&num, &DmgIid) < 0) {
        opserr << "WARNING invalid damage model ID for spring I\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (DmgIid != 0 && MatI != 0) {
        DmgI = OPS_getDamageModel(DmgIid);

        if (DmgI == NULL) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgIid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        DmgI = NULL;

      int DmgJid;
      num = 1;
      if (OPS_GetIntInput(&num, &DmgJid) < 0) {
        opserr << "WARNING invalid damage model ID for spring J\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (DmgJid != 0 && MatJ != 0) {
        DmgJ = OPS_getDamageModel(DmgJid);

        if (DmgJ == NULL) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgJid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        DmgJ = NULL;

      int DmgKid;
      num = 1;
      if (OPS_GetIntInput(&num, &DmgKid) < 0) {
        opserr << "WARNING invalid damage model ID for spring K\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (DmgKid != 0 && MatK != 0) {
        DmgK = OPS_getDamageModel(DmgKid);

        if (DmgK == NULL) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgKid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        DmgK = NULL;

      int DmgLid;
      num = 1;
      if (OPS_GetIntInput(&num, &DmgLid) < 0) {
        opserr << "WARNING invalid damage model ID for spring L\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (DmgLid != 0 && MatL != 0) {
        DmgL = OPS_getDamageModel(DmgLid);

        if (DmgL == NULL) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        DmgL = NULL;

      int PanelDmgId;
      num = 1;
      if (OPS_GetIntInput(&num, &PanelDmgId) < 0) {
        opserr << "WARNING invalid panel DmgID\n";
        opserr << "Joint2D element: " << Joint2DId << endln;
        return 0;
      }

      if (PanelDmgId != 0 && PanelMaterial != 0) {
        PanelDamage = OPS_getDamageModel(PanelDmgId);

        if (PanelDamage == NULL) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << PanelDmgId;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return 0;
        }
      } else
        DmgL = NULL;
    }

    // Create the new material
    DamageModel *damageModels[5] = {DmgI, DmgJ, DmgK, DmgL, PanelDamage};
    UniaxialMaterial *springModels[5] = {MatI, MatJ, MatK, MatL, PanelMaterial};
    theJoint2D =
        new Joint2D(Joint2DId, iNode, jNode, kNode, lNode, CenterNodeTag,
                    springModels, theDomain, LargeDisp, damageModels);
    return theJoint2D;
  } else {
    return 0;
  }
}
