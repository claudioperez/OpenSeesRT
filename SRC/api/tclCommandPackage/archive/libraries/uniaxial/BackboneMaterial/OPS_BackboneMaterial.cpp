

#include <g3_api.h>

#include <SRC/material/uniaxial/BackboneMaterial.h>

HystereticBackbone *OPS_getHystereticBackbone(int tag);

void *OPS_Backbone(G3_Runtime *rt)
{

  int argc = OPS_GetNumRemainingInputArgs() + 2;
  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Backbone tag? bbTag?\n";
    return 0;
  }

  // tag, bbTag;
  int idata[2];
  int numdata = 2;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid tags\n";
    opserr << "Backbone material: " << idata[0] << "\n";
    return 0;
  }

  HystereticBackbone *backbone = OPS_getHystereticBackbone(idata[1]);

  if (backbone == 0) {
    opserr << "WARNING backbone does not exist\n";
    opserr << "backbone: " << idata[1];
    opserr << "\nuniaxialMaterial Backbone: " << idata[0] << "\n";
    return 0;
  }

  return new BackboneMaterial(idata[0], *backbone);
}
