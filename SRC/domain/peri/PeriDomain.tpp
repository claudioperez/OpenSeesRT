#include <cmath>


template <int ndim>
PeriDomain<ndim>::PeriDomain(int totnode, int maxfam)
  // Call the constructors for our member data
  : PeriDomainBase(totnode, maxfam), pts(totnode) // initialize the container of particles with size `totnode`
{
    // Initialize each particle and set all its data fields to zero
    for (PeriParticle<ndim>& node : pts) {
        node.numfam = 0;                 // Initialize number of families to zero
        node.nodefam.assign(maxfam, 0);  // Assuming maxfam is the size of nodefam
        node.vol.assign(maxfam, 0.0);    // Volume of every neighbor, hence size maxfam
        node.vol_h = 0.0;                // Volume of the whole horizon
        node.is_force_bound.fill(0);     // Boundary condition flag for force
        node.is_disp_bound.fill(0);      // Boundary condition flag for displacement
        node.correction.assign(maxfam, 1.0); // Correction factor for each neighbor
        node.bond_dmg.assign(maxfam, 0.0); // Initialize bond damage with maxfam size

        node.coord.fill(0.0);            // coordinates in 3D space
        node.disp.fill(0.0);             // Displacement in 3D, hence 3 components
        node.vel.fill(0.0);              // Velocity in 3D
        node.acc.fill(0.0);              // Acceleration in 3D
        node.bforce.fill(0.0);           // Boundary force for each particle, so 3 components
        node.bdisp.fill(0.0);            // Boundary displacement, similar to bond force
        node.pforce.fill(0.0);           // PD force density, 3 components for 3D
        node.stress.fill(0.0);           // Assuming stress is a tensor with 6 components (symmetric)
        node.strain.fill(0.0);           // Similarly, strain is also 6 components
        node.energy = 0.0;               // Initialize scalar energy to zero
        node.misc.fill(0.0);             // Miscellaneous field (size is 6 for now)
    }
}

template <int ndim>
void PeriDomain<ndim>::set_coord(int i, const std::array<double, ndim>& coord) {
    pts[i].coord = coord; // Set the coordinates of the particle at index i
}

template <int ndim>
void PeriDomain<ndim>::create_fam(const double delta) {
    // Set the size of horizon delta
    this->delta = delta;
    // Create families for each particle
    for (int i = 0; i < totnode; i++) {
        for (int j = i+1; j < totnode; j++) {
            // Calculate the distance between the particles
            double dist = 0.0;
            for (int k = 0; k < ndim; k++) {
                dist += (pts[j].coord[k] - pts[i].coord[k]) * (pts[j].coord[k] - pts[i].coord[k]);
            }
            dist = std::sqrt(dist);
            // If the distance is less than delta, add the particle to the family
            if (dist < delta && dist > 1.0e-8*delta) {
                pts[i].nodefam[pts[i].numfam] = j;
                pts[i].numfam++;
                pts[j].nodefam[pts[j].numfam] = i;
                pts[j].numfam++;
            }
        }
    }
}

template <int ndim>
void PeriDomain<ndim>::set_vols(int i, double vol_i) {
    int numfam = pts[i].numfam;
    // Set the volume of the particle i in the end of family list
    // because in c++, indexing starts from 0
    // so the last element of the array is at index numfam-1
    // hence, the volume of the particle i is at index numfam
    pts[i].vol[numfam] = vol_i;
}

template <int ndim>
void PeriDomain<ndim>::calc_vols(double space) {
    
    int j;
    double fam_vol_i = 0.0, dist = 0.0;
    double vol_init = 0.0, coeff_mod = 0.0;

    this->space = space;

    // calculate the volume of horizon of particle i
    for (int i = 0; i < totnode; i++) {

        fam_vol_i = pts[i].vol[pts[i].numfam];
        
        for (int ind = 0; ind < pts[i].numfam; ind++) {
            j = pts[i].nodefam[ind];
            vol_init = pts[j].vol[pts[j].numfam];
            // calculate the distance between the particles i and j
            dist = 0.0;
            for (int k = 0; k < ndim; k++) {
                dist += (pts[j].coord[k] - pts[i].coord[k]) * (pts[j].coord[k] - pts[i].coord[k]);
            }
            dist = std::sqrt(dist);
            // modify the volume of the particle j in particle i's horizon
            if (dist < (this->delta - 0.5*this->space)) {
                coeff_mod = 1.0;
            }else if (dist < (this->delta + 0.5*this->space)) {
                coeff_mod = (this->delta + 0.5*this->space - dist)/this->space;
            }else {
                coeff_mod = 0.0;
            }
            pts[i].vol[ind] = vol_init * coeff_mod;
            fam_vol_i += pts[i].vol[ind];
        }

        pts[i].vol_h = fam_vol_i;
    }
}

template <int ndim>
void PeriDomain<ndim>::calc_surf_correction(){

    int j;
    double vol_h_max = 0.0, tmp_corr = 1.0;

    for (int i = 0; i < totnode; i++) {
        if (pts[i].vol_h > vol_h_max) {
            vol_h_max = pts[i].vol_h;
        }
    }

    for (int i = 0; i < totnode; i++) {
        for (int ind = 0; ind < pts[i].numfam; ind++) {
            j = pts[i].nodefam[ind];
            if (i < j) {
                // set stiffness correction factor for this bond
                tmp_corr = 2.0 * vol_h_max / (pts[i].vol_h + pts[j].vol_h);
                pts[i].correction[ind] *= tmp_corr;
                // also set for j, since it is symmetric
                for (int ind2 = 0; ind2 < pts[j].numfam; ind2++) {
                    if (pts[j].nodefam[ind2] == i) {
                        pts[j].correction[ind2] *= tmp_corr;
                        break;
                    }
                }
            }
        }
    }
}


template <int ndim>
void PeriDomain<ndim>::break_bond(const int node1, const int node2) {
    
    int nodej;
    // Remove this neighbor by replacing it with the last neighbor
    // on the list, remove the last neighbor on the list by
    // replacing it with -1, then reducing the number of neighbors by 1
    for (int i = 0; i < pts[node1].numfam; i++) {
        nodej = pts[node1].nodefam[i];
        if (nodej == node2) {
            pts[node1].nodefam[i] = pts[node1].nodefam[pts[node1].numfam-1];
            pts[node1].nodefam[pts[node1].numfam-1] = -1;
            pts[node1].numfam--;
            break;
        }
    }
    // do the same thing for the other node
    for (int i = 0; i < pts[node2].numfam; i++) {
        nodej = pts[node2].nodefam[i];
        if (nodej == node1) {
            pts[node2].nodefam[i] = pts[node2].nodefam[pts[node2].numfam-1];
            pts[node2].nodefam[pts[node2].numfam-1] = -1;
            pts[node2].numfam--;
            break;
        }
    }
}
