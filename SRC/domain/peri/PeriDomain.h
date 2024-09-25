#include <vector>
#include <array>
#include <PeriParticle.h>
#include <PeriDomainBase.h>

// ============================================
// VectorND and MatrixND are not required here
// because the linear algebra operations are not used
// --------------------------------------------
// #include <Matrix.h>
// #include <Vector.h>
// #include <VectorND.h>
// #include <MatrixND.h>
// using OpenSees::VectorND;
// using OpenSees::MatrixND;
// ============================================

template <int ndim>
class PeriDomain : public PeriDomainBase {
public:
  // ============================================
  // MEMBER FUNCTIONS
  // ============================================

  // the special "constructor" member function; no return type declared
  PeriDomain(int totnode, int maxfam);
  // override the the getNDM function in PeriDomainBase
  virtual int getNDM() const { return ndim; } // Return the number of dimensions

  void set_coord(int i, const std::array<double, ndim>& coord); // Set the coordinates of the particle at index i

  // ============================================
  // MEMBER DATA
  // ============================================
  int       totnode, maxfam;
  char      plane_type;
  std::vector<PeriParticle<ndim>> pts; // container of particles
};

#include <PeriDomain.tpp> // Include the implementation of the template class