//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#ifndef FrameTransform_h
#define FrameTransform_h

#include <vector>
#include <stdexcept>
#include <VectorND.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <CrdTransf.h>

#define MAYBE_STATIC static

using OpenSees::VectorND;
using OpenSees::MatrixND;
using OpenSees::Matrix3D;

typedef std::vector<int> Layout;

enum {
 CRDTR_TAG_CorotFrameTransfWarping3d,
 CRDTR_TAG_CorotFrameTransf3d,
 CRDTR_TAG_LinearFrameTransf3d,
 CRDTR_TAG_PDeltaFrameTransf3d
};

//
// 2D
//
class FrameTransform2d : public CrdTransf {
  public:
  FrameTransform2d(int tag, int classTag) : CrdTransf(tag, classTag) {};
};

//
// 3D
//
class FrameTransform3d : public CrdTransf {
public:
  enum {
    N, Vy, Vz, T, My, Mz,
  };

public:
  FrameTransform3d(int tag, int classTag) : CrdTransf(tag, classTag) {}


  // TODO(cmp) : make almost everything pure virtual
  virtual FrameTransform3d *getCopy() {
    return nullptr;
  }

  virtual CrdTransf *getCopy3d() {
    return getCopy();
  }


  /*
  virtual VectorND<ndf>  pullNodeUnknowns(int tag, int rate);
  virtual VectorND<ndm>  pullNodePosition(int tag, State state);
  virtual Rotation       pullNodeRotation(int tag, State state);
  virtual VectorND<ndm>  pullNodeVelocity(int tag);
  virtual VectorND<ndm>  pullNodeLocation(int tag, State state);
  */

//template <int n>
//VectorND<n> pushResponse(const VectorND<n>& response) {
//  VectorND<n> pushed;
//  const Vector wrap(response);
//  pushResponse(Vector(pushed), response);
//  return pushed;
//}

  //
  //
  //

  virtual VectorND<12>    pushResponse(VectorND<12>&pl) {
    static VectorND<12> empty{};
    return empty;
  }

  virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const {
    static VectorND<12> empty{};
    return empty;
  }
  virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) {
    static MatrixND<12,12> empty{};
    return empty;
  }
  virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) {
    static MatrixND<12,12> empty{};
    return empty;
  }

#if 0
  template<int nn, int ndf=6>
  VectorND<nn*ndf> push(const VectorND<nn*ndf>&q) {
      constexpr int n = nn*ndf;

      constexpr static int layout[] {
        N, Vy, Vz, T, My, Mz
      };

      VectorND<n> p;

      const Layout &u_layout = this->getForceLayout();
      const Layout &n_layout = this->getNodeLayout();
      int m = u_layout.size();
      Vector qt(m);
      static Vector p0(m);

      for (int node=0; node<nn; node++)
        for (int i=0; i<ndf; i++)
          for (int ti =0; ti < m; ti++)
            if (layout[i] == u_layout[ti] && node == n_layout[ti]) {
              qt[ti] = q[i];
            }
      
      const Vector& pt = this->getGlobalResistingForce(qt, p0);

      for (int node=0; node<nn; node++)
        for (int i=0; i<ndf; i++)
          for (int ti =0; ti < m; ti++)
            if (layout[i] == u_layout[ti] && node == n_layout[ti])
              p[i] = pt[ti];

      return p;
  }

  template<int nn, int ndf=6>
  MatrixND<nn*ndf,nn*ndf> push(const MatrixND<nn*ndf,nn*ndf>&k, const VectorND<nn*ndf>*q) {
    constexpr static int layout[] {
      N, Vy, Vz, T, My, Mz
    };
    constexpr int n = nn*ndf;
    MatrixND<n,n> M;

    const Layout &u_layout = this->getForceLayout();
    const Layout &n_layout = this->getNodeLayout();
    int m = u_layout.size();

    Matrix kt(m,m);

    // (1/3)
    for (int nodei = 0; nodei<nn; nodei++) {
      for (int i=0; i<ndf; i++) {
        for (int ti=0; ti<m; ti++) {
          if (layout[i] == u_layout[ti]) {
            for (int nodej = 0; nodej<nn; nodej++) {
              for (int j=0; j<ndf; j++) {
                kt(i,j) = 0.0;
                for (int tj=0; tj<m; tj++) {
                  if (layout[j] == u_layout[tj]) {
                    kt(tj,ti) = k(nodej*ndf+j,nodei*ndf+i);
                  }
                }
              }
            }
          }
        }
      }
    }

    // (2/3) Perform the computation
    const Matrix& Mt = (q == nullptr)
                       ? this->getGlobalMatrixFromLocal(kt)
                       : this->getGlobalStiffMatrix(kt, q);

    // (3/3) Transform back
    for (int nodei = 0; nodei<nn; nodei++) {
      for (int i=0; i<ndf; i++) {
        for (int ti=0; ti<m; ti++) {
          if (layout[i] == u_layout[ti] && nodei == n_layout[ti]) {
            for (int nodej = 0; nodej<nn; nodej++) {
              for (int j=0; j<n; j++) {
                M(i,j) = 0.0;
                for (int tj=0; tj<m; tj++) {
                  if (layout[j] == u_layout[tj] && nodej == n_layout[tj]) {
                    M(j+nodej*ndf,i+nodei*ndf) = Mt(tj,ti);
                  }
                }
              }
            }
          }
        }
      }
    }
    return M;
  }
#endif

protected:
  virtual const Layout& getForceLayout() const {
    static std::vector<int> l(0);
    return l;
  }

  virtual const Layout& getNodeLayout() const {
    static std::vector<int> l(0);
    return l;
  }

};

static inline VectorND<12>
getLocal(double ug[12], const Matrix3D& R, double nodeIOffset[], double nodeJOffset[])
{
  VectorND<12> ul;

  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      ul[i*3+j] = R(0,j)*ug[3*i] + R(1,j)*ug[3*i+1] + R(2,j)*ug[3*i+2];

  double Wu[3];
  if (nodeIOffset) {
    Wu[0] =  nodeIOffset[2] * ug[4] - nodeIOffset[1] * ug[5];
    Wu[1] = -nodeIOffset[2] * ug[3] + nodeIOffset[0] * ug[5];
    Wu[2] =  nodeIOffset[1] * ug[3] - nodeIOffset[0] * ug[4];

    ul[0] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[1] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[2] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  if (nodeJOffset) {
    Wu[0] =  nodeJOffset[2] * ug[10] - nodeJOffset[1] * ug[11];
    Wu[1] = -nodeJOffset[2] * ug[ 9] + nodeJOffset[0] * ug[11];
    Wu[2] =  nodeJOffset[1] * ug[ 9] - nodeJOffset[0] * ug[10];

    ul[6] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[7] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[8] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  return ul;
}

static inline VectorND<6>
getBasic(double ug[12], const Matrix3D& R, double nodeIOffset[], double nodeJOffset[], double oneOverL)
{
  VectorND<6> ub;
  VectorND<12> ul = getLocal(ug, R, nodeIOffset, nodeJOffset);

#if 0
  static double ul[12];

  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      ul[i*3+j] = R(0,j)*ug[3*i] + R(1,j)*ug[3*i+1] + R(2,j)*ug[3*i+2];

  double Wu[3];
  if (nodeIOffset) {
    Wu[0] =  nodeIOffset[2] * ug[4] - nodeIOffset[1] * ug[5];
    Wu[1] = -nodeIOffset[2] * ug[3] + nodeIOffset[0] * ug[5];
    Wu[2] =  nodeIOffset[1] * ug[3] - nodeIOffset[0] * ug[4];

    ul[0] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[1] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[2] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  if (nodeJOffset) {
    Wu[0] =  nodeJOffset[2] * ug[10] - nodeJOffset[1] * ug[11];
    Wu[1] = -nodeJOffset[2] * ug[ 9] + nodeJOffset[0] * ug[11];
    Wu[2] =  nodeJOffset[1] * ug[ 9] - nodeJOffset[0] * ug[10];

    ul[6] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[7] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[8] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }
#endif

  ub[0] = ul[6] - ul[0];

  double tmp;
  tmp   = oneOverL * (ul[1] - ul[7]);
  ub[1] = ul[ 5] + tmp;
  ub[2] = ul[11] + tmp;
  tmp   = oneOverL * (ul[8] - ul[2]);
  ub[3] = ul[ 4] + tmp;
  ub[4] = ul[10] + tmp;
  ub[5] = ul[ 9] - ul[3];

  return ub;
}

static VectorND<12>
pushLocal(const Vector& q, double L)
{
  
  VectorND<12> pl;

  double q0 = q(0);
  double q1 = q(1);
  double q2 = q(2);
  double q3 = q(3);
  double q4 = q(4);
  double q5 = q(5);

  double oneOverL = 1.0 / L;

  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  = q3;
  pl[5]  = q1;
  pl[6]  = q0;                     // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  = q5;                     // Tj
  pl[10] = q4;
  pl[11] = q2;

  return pl;
}


static inline void 
formOffsets(const Matrix3D& R, 
            const double nodeIOffset[3], 
            const double nodeJOffset[3], 
            double RWI[3][3], double RWJ[3][3])
{
  if (nodeIOffset) {
    // Compute RWI
    RWI[0][0] = -R(1,0) * nodeIOffset[2] + R(2,0) * nodeIOffset[1];
    RWI[1][0] = -R(1,1) * nodeIOffset[2] + R(2,1) * nodeIOffset[1];
    RWI[2][0] = -R(1,2) * nodeIOffset[2] + R(2,2) * nodeIOffset[1];

    RWI[0][1] =  R(0,0) * nodeIOffset[2] - R(2,0) * nodeIOffset[0];
    RWI[1][1] =  R(0,1) * nodeIOffset[2] - R(2,1) * nodeIOffset[0];
    RWI[2][1] =  R(0,2) * nodeIOffset[2] - R(2,2) * nodeIOffset[0];

    RWI[0][2] = -R(0,0) * nodeIOffset[1] + R(1,0) * nodeIOffset[0];
    RWI[1][2] = -R(0,1) * nodeIOffset[1] + R(1,1) * nodeIOffset[0];
    RWI[2][2] = -R(0,2) * nodeIOffset[1] + R(1,2) * nodeIOffset[0];
  }

  if (nodeJOffset) {
    // Compute RWJ
    RWJ[0][0] = -R(1,0) * nodeJOffset[2] + R(2,0) * nodeJOffset[1];
    RWJ[1][0] = -R(1,1) * nodeJOffset[2] + R(2,1) * nodeJOffset[1];
    RWJ[2][0] = -R(1,2) * nodeJOffset[2] + R(2,2) * nodeJOffset[1];

    RWJ[0][1] =  R(0,0) * nodeJOffset[2] - R(2,0) * nodeJOffset[0];
    RWJ[1][1] =  R(0,1) * nodeJOffset[2] - R(2,1) * nodeJOffset[0];
    RWJ[2][1] =  R(0,2) * nodeJOffset[2] - R(2,2) * nodeJOffset[0];

    RWJ[0][2] = -R(0,0) * nodeJOffset[1] + R(1,0) * nodeJOffset[0];
    RWJ[1][2] = -R(0,1) * nodeJOffset[1] + R(1,1) * nodeJOffset[0];
    RWJ[2][2] = -R(0,2) * nodeJOffset[1] + R(1,2) * nodeJOffset[0];
  }
}



static MatrixND<3,3>
FrameOrientationGradient(const Vector3D& xi, const Vector3D& xj, const Vector3D& vz, int di, int dj, int dv)
{
    Vector3D v1  = xj - xi;
    double L     = v1.norm();
    Vector3D e1  = v1/L;

    Vector3D v2  = vz.cross(e1);

    Vector3D e2 = v2 / v2.norm();
//  Vector3D v3 = e1.cross(e2);
//  Vector3D e3 = v3 / v3.norm();

    //
    Vector3D dvz{0.0};
    Vector3D dxi{0.0};
    Vector3D dxj{0.0};

    if (di != 0)
      dxi(di-1) = 1.0;
    if (dj != 0)
      dxj(dj-1) = 1.0;
    if (dv != 0)
      dvz(dv-1) = 1.0;

    //

    double   dL  = 1/L*(xj - xi).dot(dxj - dxi);
    Vector3D dv1 = dxj - dxi;
    Vector3D de1 = 1/(L*L)*(dv1*L - v1*dL);

    double L2    = v2.norm();
    Vector3D dv2 = dvz.cross(e1) + vz.cross(de1);
    double dL2   = 1/L2*v2.dot(dv2);
    Vector3D de2 = 1/(L2*L2)*(dv2*L2 - v2*dL2);

    Vector3D de3 = de1.cross(e2) + e1.cross(de2);

    MatrixND<3,3> dR;
    dR(0,0) = de1(0);
    dR(1,0) = de1(1);
    dR(2,0) = de1(2);

    dR(0,1) = de2(0);
    dR(1,1) = de2(1);
    dR(2,1) = de2(2);

    dR(0,2) = de3(0);
    dR(1,2) = de3(1);
    dR(2,2) = de3(2);

    return dR;

//  return np.stack([de1,de2,de3])
}

#endif // include guard
