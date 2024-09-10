#include <stdio.h>
#include <PeriDomain.h>
#include <Information.h>
  

PeriDomain::PeriDomain(int ndim, int totnode, int maxfam)
  // Call the constructors for our member data
  : ndim(ndim), totnode(totnode), maxfam(maxfam), pts(totnode) // initialize the container of particles with size `totnode`
{
    // Initialize each particle and set all its data fields to zero
    for (int i = 0; i < totnode; ++i) {
        pts[i] = PeriParticle();  // Default constructor for PeriParticle

        // Initialize all vector fields in PeriParticle to zeros
        pts[i].coord.assign(3, 0.0);        // Assuming coord for 3D space
        pts[i].numfam = 0;                  // Initialize number of families to zero
        pts[i].nodefam.assign(maxfam, 0);   // Assuming maxfam is the size of nodefam
        pts[i].vol.assign(maxfam, 0.0);     // Volume of every neighbor, hence size maxfam
        pts[i].vol_h = 0.0;                 // Volume of the whole horizon
        pts[i].is_force_bound.assign(3, 0); // Boundary condition flag for force
        pts[i].is_disp_bound.assign(3, 0);  // Boundary condition flag for displacement
        pts[i].bforce.assign(3, 0.0);       // Boundary force for each particle, so 3 components
        pts[i].bdisp.assign(3, 0.0);        // Boundary displacement, similar to bond force
        pts[i].correction.assign(maxfam, 1.0); // Correction factor for each neighbor
        pts[i].bond_dmg.assign(maxfam, 0.0);// Initialize bond damage with maxfam size
        pts[i].pforce.assign(3, 0.0);       // PD force density, 3 components for 3D
        pts[i].stress.assign(6, 0.0);       // Assuming stress is a tensor with 6 components (symmetric)
        pts[i].strain.assign(6, 0.0);       // Similarly, strain is also 6 components
        pts[i].disp.assign(3, 0.0);         // Displacement in 3D, hence 3 components
        pts[i].vel.assign(3, 0.0);          // Velocity in 3D
        pts[i].acc.assign(3, 0.0);          // Acceleration in 3D
        pts[i].energy = 0.0;                // Initialize scalar energy to zero
        pts[i].misc.assign(6, 0.0);         // Miscellaneous field (size is up to you)
    }
}

PeriDomain::PeriDomain(int ndim, int totnode, int maxfam, char plane_type)
  // Call the constructors for our member data
  : ndim(ndim), totnode(totnode), maxfam(maxfam), plane_type(plane_type), pts(totnode) 
  // initialize the container of particles with size `totnode`
{
    // Initialize each particle and set all its data fields to zero
    for (int i = 0; i < totnode; ++i) {
        pts[i] = PeriParticle();  // Default constructor for PeriParticle

        // Initialize all vector fields in PeriParticle to zeros
        pts[i].coord.assign(2, 0.0);        // Assuming coord for 2D space
        pts[i].numfam = 0;                  // Initialize number of families to zero
        pts[i].nodefam.assign(maxfam, 0);   // Assuming maxfam is the size of nodefam
        pts[i].vol.assign(maxfam, 0.0);     // Volume of every neighbor, hence size maxfam
        pts[i].vol_h = 0.0;                 // Volume of the whole horizon
        pts[i].is_force_bound.assign(2, 0); // Boundary condition flag for force
        pts[i].is_disp_bound.assign(2, 0);  // Boundary condition flag for displacement
        pts[i].bforce.assign(2, 0.0);       // Boundary force for each particle, so 2 components
        pts[i].bdisp.assign(2, 0.0);        // Boundary displacement, similar to bond force
        pts[i].correction.assign(maxfam, 1.0); // Correction factor for each neighbor
        pts[i].bond_dmg.assign(maxfam, 0.0);// Initialize bond damage with maxfam size
        pts[i].pforce.assign(2, 0.0);       // PD force density, 2 components for 2D
        pts[i].stress.assign(3, 0.0);       // Assuming stress is a tensor with 3 components (symmetric)
        pts[i].strain.assign(3, 0.0);       // Similarly, strain is also 3 components
        pts[i].disp.assign(2, 0.0);         // Displacement in 2D, hence 2 components
        pts[i].vel.assign(2, 0.0);          // Velocity in 2D
        pts[i].acc.assign(2, 0.0);          // Acceleration in 2D
        pts[i].energy = 0.0;                // Initialize scalar energy to zero
        pts[i].misc.assign(6, 0.0);         // Miscellaneous field (size is up to you)
    }
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
