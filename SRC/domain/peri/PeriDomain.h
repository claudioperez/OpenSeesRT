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
class PeriDomain : public PeriDomainBase
{
public:
    // ============================================
    // MEMBER FUNCTIONS
    // ============================================

    // the special "constructor" member function; no return type declared
    PeriDomain(int totnode, int maxfam);
    // override the the getNDM function in PeriDomainBase
    virtual int getNDM() const { return ndim; } // Return the number of dimensions

    void set_coord(int i, const std::array<double, ndim> &coord); // Set the coordinates of the particle at index i

    void create_fam(const double delta); // Create families for each particle

    void set_vols(int i, double vol_i); // Set the volume of the particle i

    void calc_vols(const double space); // Calculate the volume of the horizons

    void calc_surf_correction(); // Calculate the surface correction

    void break_bond(const int node1, const int node2); // Break the bond between two particles

    // ============================================
    // MEMBER DATA
    // ============================================
    // int       totnode, maxfam;
    // char      plane_type;
    std::vector<PeriParticle<ndim>> pts; // container of particles
};

#include <PeriDomain.tpp> // Include the implementation of the template class