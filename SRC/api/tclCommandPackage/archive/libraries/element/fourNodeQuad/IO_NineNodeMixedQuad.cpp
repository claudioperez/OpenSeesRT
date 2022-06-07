

#include <g3_api.h>


#include <SRC/element/fourNodeQuad/NineNodeMixedQuad.h>
void *OPS_NineNodeMixedQuad()
{
  if (OPS_GetNDM() != 2 || OPS_GetNDF() != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return 0;
  }

  // check the number of arguments is correct
  if (OPS_GetNumRemainingInputArgs() < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element NineNodeMixedQuad  eleTag?"
           << " iNode? jNode? kNode? lNode? mNode, nNode, pNode, qNode, "
              "centerNode "
           << " matTag?\n";
    return 0;
  }

  // get the id and end nodes
  int idata[11];
  int numdata = 11;
  // int NineNodeMixedQuadId, iNode, jNode, kNode, lNode ;
  // int                  mNode, nNode, pNode, qNode ;
  // int centerNode ;
  // int matID;

  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid NineNodeMixedQuad int inputs" << endln;
    return 0;
  }

  NDMaterial *theMaterial = OPS_getNDMaterial(idata[10]);

  if (theMaterial == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << idata[10];
    opserr << "\nNineNodeMixedQuad element: " << idata[0] << endln;
    return 0;
  }

  // now create the NineNodeMixedQuad and add it to the Domain
  NineNodeMixedQuad *theNineNodeMixed = new NineNodeMixedQuad(
      idata[0], idata[1], idata[2], idata[3], idata[4], idata[5], idata[6],
      idata[7], idata[8], idata[9], *theMaterial);

  return theNineNodeMixed;
}
