

#include <g3_api.h>


#include <SRC/element/updatedLagrangianBeamColumn/Inelastic2DYS03.h>
void *OPS_Inelastic2DYS03()
{
  if (OPS_GetNumRemainingInputArgs() < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dYS03 tag? Nd1? Nd2? A_ten? A_com? E? IzPos? "
              "IzNeg? ysID1? ysID2? algo?";

    return 0;
  }

  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid element2dYS int inputs" << endln;
    return 0;
  }

  int tag = idata[0];
  int ndI = idata[1];
  int ndJ = idata[2];

  double data[5];
  numdata = 5;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid element2dYS double inputs" << endln;
    return 0;
  }

  double aTens = data[0];
  double aComp = data[1];
  double E = data[2];
  double Ipos = data[3];
  double Ineg = data[4];

  numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid element2dYS int inputs" << endln;
    return 0;
  }
  int ysID1 = idata[0];
  int ysID2 = idata[1];
  int rf_algo = idata[2];

  YieldSurface_BC *theYS1 = OPS_getYieldSurface_BC(ysID1);
  if (theYS1 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID1 << endln;
    return 0;
  }

  YieldSurface_BC *theYS2 = OPS_getYieldSurface_BC(ysID2);
  if (theYS2 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID2 << endln;
    return 0;
  }

  return new Inelastic2DYS03(tag, aTens, aComp, E, Ipos, Ineg, ndI, ndJ, theYS1,
                             theYS2, rf_algo);
}
