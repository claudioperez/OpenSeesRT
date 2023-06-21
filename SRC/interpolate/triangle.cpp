double
Tri31::shapeFunction(double xi, double eta)
{
  const Vector &nd1Crds = theNodes[0]->getCrds();
  const Vector &nd2Crds = theNodes[1]->getCrds();
  const Vector &nd3Crds = theNodes[2]->getCrds();

  shp[2][0] = xi;           // N_1
  shp[2][1] = eta;          // N_2
  shp[2][2] = 1 - xi - eta; // N_3

  double J[2][2];

  // See p 180 "A First Course in Finite Elements" by Fish and Belytschko.
  J[0][0] = (nd1Crds(0) - nd3Crds(0));
  J[0][1] = (nd2Crds(0) - nd3Crds(0));
  J[1][0] = (nd1Crds(1) - nd3Crds(1));
  J[1][1] = (nd2Crds(1) - nd3Crds(1));

  double detJ        = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  double oneOverdetJ = 1.0 / detJ;
  double L[2][2];

  // L = inv(J)
  L[0][0] =  J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] =  J[0][0] * oneOverdetJ;

  // See Cook, Malkus, Plesha p. 169 for the derivation of these terms
  shp[0][0] = L[0][0];              // N_1,1
  shp[0][1] = L[0][1];              // N_2,1
  shp[0][2] = -(L[0][0] + L[0][1]); // N_3,1

  shp[1][0] = L[1][0];              // N_1,2
  shp[1][1] = L[1][1];              // N_2,2
  shp[1][2] = -(L[1][0] + L[1][1]); // N_3,2

  return detJ;
}

double
SixNodeTri::shapeFunction(double s, double t)
{
  const Vector &nd1Crds = theNodes[0]->getCrds();
  const Vector &nd2Crds = theNodes[1]->getCrds();
  const Vector &nd3Crds = theNodes[2]->getCrds();
  const Vector &nd4Crds = theNodes[3]->getCrds();
  const Vector &nd5Crds = theNodes[4]->getCrds();
  const Vector &nd6Crds = theNodes[5]->getCrds();

  shp[2][0] = s * (2 * s - 1);
  shp[2][1] = t * (2 * t - 1);
  shp[2][2] = (1 - s - t) * (1 - 2 * s - 2 * t);
// shp[2][2] = 1 - 3*s - 3*t + 2*s*s + 2*t*t + 4*s*t;
  shp[2][3] = 4 * s * t;
  shp[2][4] = 4 * t * (1 - s - t);
  shp[2][5] = 4 * s * (1 - s - t);

  // derivatives
  double N11 = 4 * s - 1;
  double N12 = 0;
  double N21 = 0;
  double N22 = 4 * t - 1;
  double N31 = -3 + 4 * s + 4 * t;
  double N32 = -3 + 4 * t + 4 * s;
  double N41 = 4 * t;
  double N42 = 4 * s;
  double N51 = -4 * t;
  double N52 = 4 - 4 * s - 8 * t;
  double N61 = 4 - 4 * t - 8 * s;
  double N62 = -4 * s;

  double J[2][2];

  J[0][0] = nd1Crds(0) * N11 + nd2Crds(0) * N21 + nd3Crds(0) * N31
          + nd4Crds(0) * N41 + nd5Crds(0) * N51 + nd6Crds(0) * N61;
  J[0][1] = nd1Crds(0) * N12 + nd2Crds(0) * N22 + nd3Crds(0) * N32
          + nd4Crds(0) * N42 + nd5Crds(0) * N52 + nd6Crds(0) * N62;
  J[1][0] = nd1Crds(1) * N11 + nd2Crds(1) * N21 + nd3Crds(1) * N31
          + nd4Crds(1) * N41 + nd5Crds(1) * N51 + nd6Crds(1) * N61;
  J[1][1] = nd1Crds(1) * N12 + nd2Crds(1) * N22 + nd3Crds(1) * N32
          + nd4Crds(1) * N42 + nd5Crds(1) * N52 + nd6Crds(1) * N62;

  double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

  double oneOverdetJ = 1 / detJ;
  double L[2][2];

  // L = inv(J)
  L[0][0] = J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] = J[0][0] * oneOverdetJ;

  double L00 = L[0][0];
  double L10 = L[1][0];
  double L01 = L[0][1];
  double L11 = L[1][1];

  shp[0][0] = L00 * N11 + L01 * N12;
  shp[0][1] = L00 * N21 + L01 * N22;
  shp[0][2] = L00 * N31 + L01 * N32;
  shp[0][3] = L00 * N41 + L01 * N42;
  shp[0][4] = L00 * N51 + L01 * N52;
  shp[0][5] = L00 * N61 + L01 * N62;

  shp[1][0] = L10 * N11 + L11 * N12;
  shp[1][1] = L10 * N21 + L11 * N22;
  shp[1][2] = L10 * N31 + L11 * N32;
  shp[1][3] = L10 * N41 + L11 * N42;
  shp[1][4] = L10 * N51 + L11 * N52;
  shp[1][5] = L10 * N61 + L11 * N62;

  return detJ;
}

// anonymous namespace for utilities
// from ASDEmbeddedNodeElement
namespace {

double
det2(const Matrix &J)
{
  return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
}

double
det3(const Matrix &J)
{
  return J(0, 0) * J(1, 1) * J(2, 2) - J(0, 0) * J(1, 2) * J(2, 1) -
         J(0, 1) * J(1, 0) * J(2, 2) + J(0, 1) * J(1, 2) * J(2, 0) +
         J(0, 2) * J(1, 0) * J(2, 1) - J(0, 2) * J(1, 1) * J(2, 0);
}

void
cross(const Vector &a, const Vector &b, Vector &c)
{
  c(0) = a(1) * b(2) - a(2) * b(1);
  c(1) = a(2) * b(0) - a(0) * b(2);
  c(2) = a(0) * b(1) - a(1) * b(0);
}

namespace tri {

double
shapeFun(double x, double y, int i)
{
  if (i == 0)
    return 1.0 - x - y;
  else if (i == 1)
    return x;
  else if (i == 2)
    return y;
  return 0.0;
}

void
shapeFunDer(Matrix &dN)
{
  dN(0, 0) = -1.0;
  dN(0, 1) = -1.0;
  dN(1, 0) = 1.0;
  dN(1, 1) = 0.0;
  dN(2, 0) = 0.0;
  dN(2, 1) = 1.0;
}

void
globalCoord(const Matrix &X, double lx, double ly, double &gx, double &gy)
{
  gx = gy = 0.0;
  for (int i = 0; i < 3; i++) {
    double N = shapeFun(lx, ly, i);
    gx += N * X(0, i);
    gy += N * X(1, i);
  }
}

void
globalCoord(const Matrix &X, double lx, double ly, double &gx, double &gy,
            double &gz)
{
  gx = gy = gz = 0.0;
  for (int i = 0; i < 3; i++) {
    double N = shapeFun(lx, ly, i);
    gx += N * X(0, i);
    gy += N * X(1, i);
    gz += N * X(2, i);
  }
}

void
localCoord(const Matrix &X, const Matrix &invJ, double gx, double gy,
           double &lx, double &ly)
{
  lx = ly = 0.0;
  double px, py;
  globalCoord(X, lx, ly, px, py);
  Vector D(2);
  Vector DL(2);
  D(0) = gx - px;
  D(1) = gy - py;
  DL.addMatrixVector(0.0, invJ, D, 1.0);
  lx = DL(0);
  ly = DL(1);
}

void
localCoord(const Matrix &X, const Matrix &invJ, double gx, double gy, double gz,
           double &lx, double &ly)
{
  lx = ly = 0.0;
  double px, py, pz;
  globalCoord(X, lx, ly, px, py, pz);
  Vector D(3);
  Vector DL(3);
  D(0) = gx - px;
  D(1) = gy - py;
  D(2) = gz - pz;
  DL.addMatrixVector(0.0, invJ, D, 1.0);
  lx = DL(0);
  ly = DL(1);
}

void
fillVzInJacobian(Matrix &J)
{
  double nx   = J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0);
  double ny   = J(0, 1) * J(2, 0) - J(0, 0) * J(2, 1);
  double nz   = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
  double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (norm > std::numeric_limits<double>::epsilon()) {
    J(0, 2) = nx / norm;
    J(1, 2) = ny / norm;
    J(2, 2) = nz / norm;
  }
}

} // namespace tri

namespace tet {

double
shapeFun(double x, double y, double z, int i)
{
  if (i == 0)
    return 1.0 - (x + y + z);
  else if (i == 1)
    return x;
  else if (i == 2)
    return y;
  else if (i == 3)
    return z;
  return 0.0;
}

void
shapeFunDer(Matrix &dN)
{
  dN(0, 0) = -1.0;
  dN(0, 1) = -1.0;
  dN(0, 2) = -1.0;
  dN(1, 0) =  1.0;
  dN(1, 1) =  0.0;
  dN(1, 2) =  0.0;
  dN(2, 0) =  0.0;
  dN(2, 1) =  1.0;
  dN(2, 2) =  0.0;
  dN(3, 0) =  0.0;
  dN(3, 1) =  0.0;
  dN(3, 2) =  1.0;
}

void
globalCoord(const Matrix &X, double lx, double ly, double lz, double &gx,
            double &gy, double &gz)
{
  gx = gy = gz = 0.0;
  for (int i = 0; i < 4; i++) {
    double N = shapeFun(lx, ly, lz, i);
    gx += N * X(0, i);
    gy += N * X(1, i);
    gz += N * X(2, i);
  }
}

void
localCoord(const Matrix &X, const Matrix &invJ, double gx, double gy, double gz,
           double &lx, double &ly, double &lz)
{
  lx = ly = lz = 0.0;
  double px, py, pz;
  globalCoord(X, lx, ly, lz, px, py, pz);
  Vector D(3);
  Vector DL(3);
  D(0) = gx - px;
  D(1) = gy - py;
  D(2) = gz - pz;
  DL.addMatrixVector(0.0, invJ, D, 1.0);
  lx = DL(0);
  ly = DL(1);
  lz = DL(2);
}

} // namespace tet

} // namespace
