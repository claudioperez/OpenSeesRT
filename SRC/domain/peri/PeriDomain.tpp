
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
}

template <int ndim>
void PeriDomain<ndim>::set_coord(int i, const std::array<double, ndim>& coord) {
    pts[i].coord = coord; // Set the coordinates of the particle at index i
}