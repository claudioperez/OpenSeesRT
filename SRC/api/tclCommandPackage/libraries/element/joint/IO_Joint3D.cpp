

#include <g3_api.h>


#include <SRC/element/joint/Joint3D.h>
void *OPS_Joint3D()
{
  if (OPS_GetNDM() != 3 || OPS_GetNDF() != 6) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with Joint3D element\n";
    return 0;
  }

  // check the number of arguments is correct
  if (OPS_GetNumRemainingInputArgs() != 12 &&
      OPS_GetNumRemainingInputArgs() != 16) {
    opserr << "WARNING incorrect number of arguments\n";
    // printCommand(argc, argv);
    opserr << "Want:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? "
              "MatX? MatY? MatZ? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? "
              "MatX? MatY? MatZ? LrgDsp? -damage DmgX DmgY DmgZ\n";
    return 0;
  }

  // get the id and end nodes
  int idata[8];
  int num = 8;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING invalid Joint3D int inputs" << endln;
    return 0;
  }
  int Joint3DId = idata[0];
  int iNode = idata[1];
  int jNode = idata[2];
  int kNode = idata[3];
  int lNode = idata[4];
  int mNode = idata[5];
  int nNode = idata[6];
  ;

  // Get the center node
  int CenterNodeTag = idata[7];

  // check domain for existence of internal node tag
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0)
    return 0;
  Node *CenterNode = theDomain->getNode(CenterNodeTag);
  if (CenterNode != 0) {
    opserr
        << "WARNING node tag specified for the center node already exists.\n";
    opserr << "Use a new node tag.\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return 0;
  }

  UniaxialMaterial *MatX = NULL;
  int MatXid;
  num = 1;
  if (OPS_GetIntInput(&num, &MatXid) < 0) {
    opserr << "WARNING invalid material ID for spring X\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return 0;
  }

  MatX = OPS_getUniaxialMaterial(MatXid);
  if (MatX == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatXid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return 0;
  }

  UniaxialMaterial *MatY = NULL;
  int MatYid;
  num = 1;
  if (OPS_GetIntInput(&num, &MatYid) < 0) {
    opserr << "WARNING invalid material ID for spring Y\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return 0;
  }

  MatY = OPS_getUniaxialMaterial(MatYid);
  if (MatY == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatYid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return 0;
  }

  UniaxialMaterial *MatZ = NULL;
  int MatZid;
  num = 1;
  if (OPS_GetIntInput(&num, &MatZid) < 0) {
    opserr << "WARNING invalid material ID for spring Z\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return 0;
  }

  MatZ = OPS_getUniaxialMaterial(MatZid);
  if (MatZ == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatZid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return 0;
  }

  int LargeDisp;
  num = 1;
  if (OPS_GetIntInput(&num, &LargeDisp) < 0) {
    // use 0 as default
    LargeDisp = 0;
  }

  Joint3D *theJoint3D;
  // Decide to use which constructor, based on the number of arguments
  if (OPS_GetNumRemainingInputArgs() == 12) {

    // Using Joint3D constructor without damage
    UniaxialMaterial *springModels[3] = {MatX, MatY, MatZ};
    theJoint3D =
        new Joint3D(Joint3DId, iNode, jNode, kNode, lNode, mNode, nNode,
                    CenterNodeTag, springModels, theDomain, LargeDisp);

    // if get here we have successfully created the element and added it to the
    // domain
    return theJoint3D;
  }

  else // if ( (argc-argStart) == 16  )
  {
    opserr << "WARNING Using Joint3D constructor with damage not implemented "
              "in this version\n";
    return 0;
  }
  return 0;
}
