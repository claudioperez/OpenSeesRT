

#include <g3_api.h>


#include <SRC/element/updatedLagrangianBeamColumn/Inelastic2DYS02.h>
void *OPS_Inelastic2DYS02()
{

  if (OPS_GetNumRemainingInputArgs() < 12) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? "
              "cycType? wt? power? algo?";

    return 0;
  }

  int idata[3];
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid element2dYS int inputs" << endln;
  }

  int tag = idata[0];
  int ndI = idata[1];
  int ndJ = idata[2];

  double data[3];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid element2dYS double inputs" << endln;
  }

  double A = data[0];
  double E = data[1];
  double I = data[2];

  numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid element2dYS int inputs" << endln;
  }

  int ysID1 = idata[0];
  int ysID2 = idata[1];
  int cyc_type = idata[2];

  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid element2dYS double inputs" << endln;
  }

  double delpmax = data[0];
  double alfa = data[1];
  double beta = data[2];

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

  int rf_algo = -1;

  CyclicModel *theModel = OPS_getCyclicModel(cyc_type);

  return new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1, theYS2, theModel,
                             delpmax, alfa, beta, rf_algo);
}
