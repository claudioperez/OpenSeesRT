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
#include <CrdTransf.h>

#define MAYBE_STATIC static

using OpenSees::VectorND;
using OpenSees::MatrixND;

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
//  throw std::runtime_error("");
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

#endif // include guard
