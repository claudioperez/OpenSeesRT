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
  PeriDomain(int totnode, int maxfam, char plane_type);



  // // A normal member function
  // int hello(double);

  // // Print a representation of the domain
  // void print(int flag);


  // ============================================
  // MEMBER DATA
  // ============================================
  int       totnode, maxfam;
  char      plane_type;
  std::vector<PeriParticle<ndim>> pts; // container of particles
};

template <int ndim>
PeriDomain<ndim>::PeriDomain(int totnode, int maxfam)
  // Call the constructors for our member data
  : PeriDomainBase(totnode, maxfam), pts(totnode) // initialize the container of particles with size `totnode`
{
    // Initialize each particle and set all its data fields to zero
    for (PeriParticle<ndim>& node : pts) {
      node.coord.fill(0.0);              // coordinates in 3D space
        node.numfam = 0;                 // Initialize number of families to zero
        node.nodefam.assign(maxfam, 0);  // Assuming maxfam is the size of nodefam
        node.vol.assign(maxfam, 0.0);    // Volume of every neighbor, hence size maxfam
        node.vol_h = 0.0;                // Volume of the whole horizon
        node.is_force_bound.fill(0);     // Boundary condition flag for force
        node.is_disp_bound.fill(0);      // Boundary condition flag for displacement
        node.bforce.fill(0.0);           // Boundary force for each particle, so 3 components
        node.bdisp.fill(0.0);            // Boundary displacement, similar to bond force
        node.correction.assign(maxfam, 1.0); // Correction factor for each neighbor
        node.bond_dmg.assign(maxfam, 0.0); // Initialize bond damage with maxfam size
        node.pforce.fill(0.0);           // PD force density, 3 components for 3D
        node.stress.fill(0.0);           // Assuming stress is a tensor with 6 components (symmetric)
        node.strain.fill(0.0);           // Similarly, strain is also 6 components
        node.disp.fill(0.0);             // Displacement in 3D, hence 3 components
        node.vel.fill(0.0);              // Velocity in 3D
        node.acc.fill(0.0);              // Acceleration in 3D
        node.energy = 0.0;               // Initialize scalar energy to zero
        node.misc.fill(0.0);             // Miscellaneous field (size is 6 for now)
    }

    // for (int i = 0; i < totnode; ++i) {
    //     pts[i] = PeriParticle<ndim>();  // Default constructor for PeriParticle

    //     // Initialize all vector fields in PeriParticle to zeros
    //     pts[i].coord.assign(3, 0.0);        // Assuming coord for 3D space
    //     pts[i].numfam = 0;                  // Initialize number of families to zero
    //     pts[i].nodefam.assign(maxfam, 0);   // Assuming maxfam is the size of nodefam
    //     pts[i].vol.assign(maxfam, 0.0);     // Volume of every neighbor, hence size maxfam
    //     pts[i].vol_h = 0.0;                 // Volume of the whole horizon
    //     pts[i].is_force_bound.assign(3, 0); // Boundary condition flag for force
    //     pts[i].is_disp_bound.assign(3, 0);  // Boundary condition flag for displacement
    //     pts[i].bforce.assign(3, 0.0);       // Boundary force for each particle, so 3 components
    //     pts[i].bdisp.assign(3, 0.0);        // Boundary displacement, similar to bond force
    //     pts[i].correction.assign(maxfam, 1.0); // Correction factor for each neighbor
    //     pts[i].bond_dmg.assign(maxfam, 0.0);// Initialize bond damage with maxfam size
    //     pts[i].pforce.assign(3, 0.0);       // PD force density, 3 components for 3D
    //     pts[i].stress.assign(6, 0.0);       // Assuming stress is a tensor with 6 components (symmetric)
    //     pts[i].strain.assign(6, 0.0);       // Similarly, strain is also 6 components
    //     pts[i].disp.assign(3, 0.0);         // Displacement in 3D, hence 3 components
    //     pts[i].vel.assign(3, 0.0);          // Velocity in 3D
    //     pts[i].acc.assign(3, 0.0);          // Acceleration in 3D
    //     pts[i].energy = 0.0;                // Initialize scalar energy to zero
    //     pts[i].misc.assign(6, 0.0);         // Miscellaneous field (size is up to you)
    // }
}

template <int ndim>
PeriDomain<ndim>::PeriDomain(int totnode, int maxfam, char plane_type)
  // Call the constructors for our member data
  :  PeriDomainBase(totnode, maxfam), plane_type(plane_type), pts(totnode) 
  // initialize the container of particles with size `totnode`
{
    // Initialize each particle and set all its data fields to zero
    // for (int i = 0; i < totnode; ++i) {
    //     pts[i] = PeriParticle<ndim>();  // Default constructor for PeriParticle

    //     // Initialize all vector fields in PeriParticle to zeros
    //     pts[i].coord.assign(2, 0.0);        // Assuming coord for 2D space
    //     pts[i].numfam = 0;                  // Initialize number of families to zero
    //     pts[i].nodefam.assign(maxfam, 0);   // Assuming maxfam is the size of nodefam
    //     pts[i].vol.assign(maxfam, 0.0);     // Volume of every neighbor, hence size maxfam
    //     pts[i].vol_h = 0.0;                 // Volume of the whole horizon
    //     pts[i].is_force_bound.assign(2, 0); // Boundary condition flag for force
    //     pts[i].is_disp_bound.assign(2, 0);  // Boundary condition flag for displacement
    //     pts[i].bforce.assign(2, 0.0);       // Boundary force for each particle, so 2 components
    //     pts[i].bdisp.assign(2, 0.0);        // Boundary displacement, similar to bond force
    //     pts[i].correction.assign(maxfam, 1.0); // Correction factor for each neighbor
    //     pts[i].bond_dmg.assign(maxfam, 0.0);// Initialize bond damage with maxfam size
    //     pts[i].pforce.assign(2, 0.0);       // PD force density, 2 components for 2D
    //     pts[i].stress.assign(3, 0.0);       // Assuming stress is a tensor with 3 components (symmetric)
    //     pts[i].strain.assign(3, 0.0);       // Similarly, strain is also 3 components
    //     pts[i].disp.assign(2, 0.0);         // Displacement in 2D, hence 2 components
    //     pts[i].vel.assign(2, 0.0);          // Velocity in 2D
    //     pts[i].acc.assign(2, 0.0);          // Acceleration in 2D
    //     pts[i].energy = 0.0;                // Initialize scalar energy to zero
    //     pts[i].misc.assign(6, 0.0);         // Miscellaneous field (size is up to you)
    // }
}


// int PeriDomain::hello(double x)
// {
//   // print the given double
//   printf("hello, %lf\n", x);

//   // return an integer
//   return 32;
// }

// void PeriDomain::print(int flag)
// {
//   printf("  \"Particles\": [\n");
//   for (PeriParticle<3>& particle : particles) {
//     printf("    ");
//     particle.print(flag);
//   }
//   printf("  ]\n");
// }
