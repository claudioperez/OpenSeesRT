/*
 * Summary of Functions
 *
 *
 * void UEL(double *RHS, double *AMATRX, double *SVARS, double *ENERGY, int NDOFEL, int NRHS, int NSVARS,
 *          double *PROPS, int NPROPS, double *COORDS, int MCRD, int NNODE, double *U, double *DU, double *V, double *A,
 *          int JTYPE, double *TIME, double DTIME, int KSTEP, int KINC, int JELEM, double *PARAMS, int NDLOAD, int *JDLTYP,
 *          double *ADLMAG, double *PREDEF, int NPREDF, int *LFLAGS, int MLVARX, double *DDLMAG, int MDLOAD, double *PNEWDT,
 *          double *JPROPS, int NJPROP, double PERIOD)
 *
 *    - Assembles the element stiffness matrix by determining whether the near-field dynamics bond is broken.
 *    - Inputs: RHS, AMATRX, SVARS, ENERGY, PROPS, COORDS, U, DU, V, A, TIME, DTIME, PARAMS, JDLTYP, ADLMAG, PREDEF, LFLAGS, DDLMAG, PNEWDT, JPROPS
 *              NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF, MLVARX, MDLOAD, NJPROP, PERIOD
 *    - Outputs: RHS, AMATRX, SVARS, ENERGY, PROPS, COORDS, U, DU, V, A, TIME, DTIME, PARAMS, JDLTYP, ADLMAG, PREDEF, LFLAGS, DDLMAG, PNEWDT, JPROPS
 *
 * void mat_zero(double *mat, int rows, int cols)
 *    - Initializes a matrix to zero.
 *    - Inputs: mat (matrix to be initialized), rows (number of rows), cols (number of columns)
 *    - Outputs: None
 *
 * void mat_mult(double *mat1, double *mat2, double *mat3, int row1, int col1, int row2, int col2)
 *    - Multiplies two matrices.
 *    - Inputs: mat1 (first matrix), mat2 (second matrix), row1 (rows in mat1), col1 (columns in mat1), row2 (rows in mat2), col2 (columns in mat2)
 *    - Outputs: mat3 (resultant matrix)
 *
 * void mat_tran(double *mat1, double *mat2, int row1, int col1)
 *    - Transposes a matrix.
 *    - Inputs: mat1 (matrix to be transposed), row1 (number of rows in mat1), col1 (number of columns in mat1)
 *    - Outputs: mat2 (transposed matrix)
 *
 * void calc_Jdet(double xRef[2][8], double i_xi, double i_eta, double j_xi, double j_eta, double *i_Jdet, double *j_Jdet)
 *    - Calculates the determinant of the Jacobi matrix.
 *    - Inputs: xRef (reference coordinates), i_xi, i_eta, j_xi, j_eta (parametric coordinates)
 *    - Outputs: i_Jdet, j_Jdet (determinants of the Jacobi matrix)
 *
 * void calc_gauss_intg(double i_xi, double i_eta, double x_ref[2][8], double j_xi, double j_eta, double block_stiff[16][16], double ele_esize, double param_1, double param_2)
 *    - Calculates the Gauss integration matrix.
 *    - Inputs: i_xi, i_eta, j_xi, j_eta (parametric coordinates), x_ref (reference coordinates), ele_esize, param_1, param_2 (parameters)
 *    - Outputs: block_stiff (Gauss integration matrix)
 *
 * void gauss_loc(double xCur[2][8], double loc_gauss[2][8])
 *    - Calculates the coordinates of Gauss points.
 *    - Inputs: xCur (current coordinates)
 *    - Outputs: loc_gauss (coordinates of Gauss points)
 */


#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_ELEMENTS 100000
#define MAX_RELATED_ELEMENTS 100

// Global variables
double ele_bond_sum_state_fem[4][1][MAX_ELEMENTS]; // FEM element damage
int    ele_bond_state_pd[1][1][MAX_RELATED_ELEMENTS][MAX_ELEMENTS]; // Related element determination
int    ele_bond_sum_state_pd[4][1][MAX_RELATED_ELEMENTS][MAX_ELEMENTS]; // Number of broken bonds in PD elements

void mat_zero(double *mat, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mat[i * cols + j] = 0.0;
        }
    }
}

void mat_mult(double *mat1, double *mat2, double *mat3, int row1, int col1, int row2, int col2) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            mat3[i * col2 + j] = 0.0;
            for (int k = 0; k < col1; k++) {
                mat3[i * col2 + j] += mat1[i * col1 + k] * mat2[k * col2 + j];
            }
        }
    }
}

void mat_tran(double *mat1, double *mat2, int row1, int col1) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col1; j++) {
            mat2[j * row1 + i] = mat1[i * col1 + j];
        }
    }
}

void calc_Jdet(double xRef[2][8], double i_xi, double i_eta, double j_xi, double j_eta, double *i_Jdet, double *j_Jdet) {
    double i_J[2][2] = {0}, 
           j_J[2][2] = {0};
    double i_N_diff[2][4], 
           j_N_diff[2][4];

    // Partial derivatives of shape functions with respect to parametric coordinates in the main element
    i_N_diff[0][0] = (i_eta - 1.0) / 4.0;
    i_N_diff[1][0] = (i_xi - 1.0) / 4.0;
    i_N_diff[0][1] = (1.0 - i_eta) / 4.0;
    i_N_diff[1][1] = -1.0 * (i_xi + 1.0) / 4.0;
    i_N_diff[0][2] = (1.0 + i_eta) / 4.0;
    i_N_diff[1][2] = (1.0 + i_xi) / 4.0;
    i_N_diff[0][3] = -1.0 * (1.0 + i_eta) / 4.0;
    i_N_diff[1][3] = (1.0 - i_xi) / 4.0;

    // Partial derivatives of shape functions with respect to parametric coordinates in the auxiliary element
    j_N_diff[0][0] = (j_eta - 1.0) / 4.0;
    j_N_diff[1][0] = (j_xi - 1.0) / 4.0;
    j_N_diff[0][1] = (1.0 - j_eta) / 4.0;
    j_N_diff[1][1] = -1.0 * (j_xi + 1.0) / 4.0;
    j_N_diff[0][2] = (1.0 + j_eta) / 4.0;
    j_N_diff[1][2] = (1.0 + j_xi) / 4.0;
    j_N_diff[0][3] = -1.0 * (1.0 + j_eta) / 4.0;
    j_N_diff[1][3] = (1.0 - j_xi) / 4.0;

    // Calculate Jacobi matrix for the main element
    for (int i = 0; i < 4; i++) {
        i_J[0][0] += i_N_diff[0][i] * xRef[0][i];
        i_J[0][1] += i_N_diff[0][i] * xRef[1][i];
        i_J[1][0] += i_N_diff[1][i] * xRef[0][i];
        i_J[1][1] += i_N_diff[1][i] * xRef[1][i];
    }

    // Calculate Jacobi matrix for the auxiliary element
    for (int i = 0; i < 4; i++) {
        j_J[0][0] += j_N_diff[0][i] * xRef[0][i + 4];
        j_J[0][1] += j_N_diff[0][i] * xRef[1][i + 4];
        j_J[1][0] += j_N_diff[1][i] * xRef[0][i + 4];
        j_J[1][1] += j_N_diff[1][i] * xRef[1][i + 4];
    }

    // Calculate the determinant of the Jacobi matrix
    *i_Jdet = i_J[0][0] * i_J[1][1] - i_J[0][1] * i_J[1][0];
    *j_Jdet = j_J[0][0] * j_J[1][1] - j_J[0][1] * j_J[1][0];
}

void calc_gauss_intg(double i_xi, double i_eta, 
                     double x_ref[2][8], double j_xi, double j_eta, 
                     double block_stiff[16][16], double ele_esize, double param_1, double param_2) 
{
    double peri_B[2][16] = {0}, 
           i_N[4], j_N[4], 
           peri_D[2][2], 
           i_realcoord[2] = {0}, 
           j_realcoord[2] = {0};

    double peri_DB[2][16], tran_peri_B[16][2];
    double i_Jdet, j_Jdet, kexi;

    calc_Jdet(x_ref, i_xi, i_eta, j_xi, j_eta, &i_Jdet, &j_Jdet);

    // Shape functions
    i_N[0] = 1.0 / 4.0 * (1.0 - i_xi) * (1.0 - i_eta);
    i_N[1] = 1.0 / 4.0 * (1.0 + i_xi) * (1.0 - i_eta);
    i_N[2] = 1.0 / 4.0 * (1.0 + i_xi) * (1.0 + i_eta);
    i_N[3] = 1.0 / 4.0 * (1.0 - i_xi) * (1.0 + i_eta);

    j_N[0] = 1.0 / 4.0 * (1.0 - j_xi) * (1.0 - j_eta) * -1.0;
    j_N[1] = 1.0 / 4.0 * (1.0 + j_xi) * (1.0 - j_eta) * -1.0;
    j_N[2] = 1.0 / 4.0 * (1.0 + j_xi) * (1.0 + j_eta) * -1.0;
    j_N[3] = 1.0 / 4.0 * (1.0 - j_xi) * (1.0 + j_eta) * -1.0;

    for (int i = 0; i < 4; i++) {
        peri_B[0][i * 2] = i_N[i];
        peri_B[0][i * 2 + 8] = j_N[i];
        peri_B[1][i * 2 + 1] = i_N[i];
        peri_B[1][i * 2 + 9] = j_N[i];
    }

    // Calculate current particle position coordinates
    for (int i = 0; i < 4; i++) {
        i_realcoord[0] += i_N[i] * x_ref[0][i];
        i_realcoord[1] += i_N[i] * x_ref[1][i];

        j_realcoord[0] += j_N[i] * x_ref[0][i + 4];
        j_realcoord[1] += j_N[i] * x_ref[1][i + 4];
    }

    // Distance between interaction points
    kexi = sqrt((i_realcoord[0] + j_realcoord[0]) * (i_realcoord[0] + j_realcoord[0]) + 
                (i_realcoord[1] + j_realcoord[1]) * (i_realcoord[1] + j_realcoord[1]));

    double L = ele_esize * param_2;
    double peri_Jdet = i_Jdet * j_Jdet * exp(-1.0 * kexi / L) * param_1;

    peri_D[0][0] = (i_realcoord[0] + j_realcoord[0]) * (i_realcoord[0] + j_realcoord[0]) * peri_Jdet;
    peri_D[0][1] = (i_realcoord[0] + j_realcoord[0]) * (i_realcoord[1] + j_realcoord[1]) * peri_Jdet;
    peri_D[1][0] = (i_realcoord[1] + j_realcoord[1]) * (i_realcoord[0] + j_realcoord[0]) * peri_Jdet;
    peri_D[1][1] = (i_realcoord[1] + j_realcoord[1]) * (i_realcoord[1] + j_realcoord[1]) * peri_Jdet;

    mat_zero((double *)peri_DB, 2, 16);
    mat_zero((double *)block_stiff, 16, 16);

    mat_mult((double *)peri_D, (double *)peri_B, (double *)peri_DB, 2, 2, 2, 16);
    mat_tran((double *)peri_B, (double *)tran_peri_B, 2, 16);
    mat_mult((double *)tran_peri_B, (double *)peri_DB, (double *)block_stiff, 16, 2, 2, 16);
}

static void
gauss_loc(double xCur[2][8], double loc_gauss[2][8]) {
    double x_Gspoint[4] = {-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626};
    double y_Gspoint[4] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
    double u_N[4];

    mat_zero((double *)loc_gauss, 2, 8);

    for (int i = 0; i < 4; i++) {
        u_N[0] = 1.0 / 4.0 * (1.0 - x_Gspoint[i]) * (1.0 - y_Gspoint[i]);
        u_N[1] = 1.0 / 4.0 * (1.0 + x_Gspoint[i]) * (1.0 - y_Gspoint[i]);
        u_N[2] = 1.0 / 4.0 * (1.0 + x_Gspoint[i]) * (1.0 + y_Gspoint[i]);
        u_N[3] = 1.0 / 4.0 * (1.0 - x_Gspoint[i]) * (1.0 + y_Gspoint[i]);

        for (int j = 0; j < 4; j++) {
            loc_gauss[0][i] += u_N[j] * xCur[0][j];
            loc_gauss[1][i] += u_N[j] * xCur[1][j];
            loc_gauss[0][i + 4] += u_N[j] * xCur[0][j + 4];
            loc_gauss[1][i + 4] += u_N[j] * xCur[1][j + 4];
        }
    }
}


#if 1
#include <math.h>
#include <stdio.h>

// void mat_zero(double matrix[][8], int rows, int cols);
// void gauss_loc(const double coords[][8], double gauss[][8]);
// void calc_gauss_intg(double x1, double y1, const double coords[][8], double x2, double y2, double block_stiff[][16], double prop1, double prop3, double prop4);
// void mat_mult(const double *A, const double *B, double *C, int rowsA, int colsA, int colsB, int rhs);

void UEL(double *RHS, double *AMATRX, double *SVARS, double *ENERGY, const int NDOFEL, const int NRHS, const int NSVARS,
         const double *PROPS, const int NPROPS, const double *COORDS, const int MCRD, const int NNODE, const double *U, const double *DU,
         const double *V, const double *A, const int JTYPE, const double *TIME, const double DTIME, const int KSTEP, const int KINC,
         const int JELEM, const double *PARAMS, const int NDLOAD, const int *JDLTYP, const double *ADLMAG, const double *PREDEF,
         const int NPREDF, const int *LFLAGS, const int MLVARX, const double *DDLMAG, const int MDLOAD, const double PNEWDT,
         const double *JPROPS, const int NJPROP, const double PERIOD) {

    const int ntens = 4;
    const double x_Gspoint[4] = {-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626};
    const double y_Gspoint[4] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};

    // Initialize global variables
    if (KINC == 1) {
        ele_bond_state_pd[1][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = 0;
        for (int i = 0; i < NSVARS; i++) {
            SVARS[i] = 0.0;
        }
        for (int i = 0; i < 4; i++) {
            ele_bond_sum_state_fem[i][1][(int)COORDS[1 * 9]] = 0.0;
            ele_bond_sum_state_pd[i][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = 0;
        }
    }

    // Initialize matrices
    double curr_gauss[2][8], ref_gauss[2][8], curr_coord[2][8], block_stiff[16][16];
    mat_zero((double *)curr_gauss, 2, 8);
    mat_zero((double *)ref_gauss, 2, 8);
    mat_zero((double *)curr_coord, 2, 8);
    mat_zero((double *)AMATRX, NDOFEL, NDOFEL);
    mat_zero((double *)block_stiff, 16, 16);
//  mat_zero(curr_gauss, 2, 8);
//  mat_zero(ref_gauss, 2, 8);
//  mat_zero(curr_coord, 2, 8);
//  mat_zero((double (*)[8])AMATRX, NDOFEL, NDOFEL);
//  mat_zero(block_stiff, 16, 16);
    double rhsmatrx[NDOFEL];
//  mat_zero((double (*)[1])rhsmatrx, NDOFEL, NRHS);
    mat_zero((double *)rhsmatrx, NDOFEL, NRHS);

    // Calculate current node coordinates
    for (int i = 0; i < 8; i++) {
        curr_coord[0][i] = COORDS[0 * NNODE + i] + U[i * 2 - 1];
        curr_coord[1][i] = COORDS[1 * NNODE + i] + U[i * 2];
    }

    // Calculate Gauss coordinates
    gauss_loc(curr_coord, curr_gauss);
    gauss_loc((const double (*)[8])COORDS, ref_gauss);

    // Determine bond breakage
    int loop_number  = 1,
        new_bre_bond = 0;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double ref_bond  = sqrt(pow(ref_gauss[0][i] - ref_gauss[0][j + 4], 2) 
                                  + pow(ref_gauss[1][i] - ref_gauss[1][j + 4], 2));
            double curr_bond = sqrt(pow(curr_gauss[0][i] - curr_gauss[0][j + 4], 2) 
                                  + pow(curr_gauss[1][i] - curr_gauss[1][j + 4], 2));

            if (SVARS[3 + loop_number] == 0.0) {
                if (ref_bond != 0.0) {
                    if (((curr_bond - ref_bond) / ref_bond) >= PROPS[2]) {
                        new_bre_bond++;
                        SVARS[3 + loop_number] = 1.0;
                    }
                }
            }
            loop_number++;
        }
    }

    // Assemble element stiffness matrix
    loop_number = 1;
    int para_build_bond = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (SVARS[3 + loop_number] != 1.0) {
                calc_gauss_intg(x_Gspoint[i], y_Gspoint[i], (const double (*)[8])COORDS, x_Gspoint[j], y_Gspoint[j], block_stiff, PROPS[1], PROPS[3], PROPS[4]);
                for (int k = 0; k < NDOFEL * NDOFEL; k++) {
                    AMATRX[k] += block_stiff[k / 16][k % 16];
                }
                para_build_bond++;
            }
            loop_number++;
        }
    }

    mat_mult(AMATRX, U, rhsmatrx, NDOFEL, NDOFEL, 1, 1);

    for (int i = 0; i < NDOFEL; i++) {
        RHS[i] = -rhsmatrx[i];
    }

    // Update state variables
    SVARS[0] = new_bre_bond;
    SVARS[19] = para_build_bond;
    SVARS[20] = 16 - para_build_bond;
    SVARS[1] = (16 - para_build_bond) / 16.0;

    // Store the number of broken bonds
    ele_bond_state_pd[1][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = 1;
    ele_bond_sum_state_pd[1][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = SVARS[20];
    ele_bond_sum_state_pd[2][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = SVARS[20];
    ele_bond_sum_state_pd[3][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = SVARS[20];
    ele_bond_sum_state_pd[4][1][(int)COORDS[2 * 9]][(int)COORDS[1 * 9]] = SVARS[20];
}


#else
void UEL(double *RHS, double *AMATRX, double *SVARS, double *ENERGY, int NDOFEL, int NRHS, int NSVARS,
         double *PROPS, int NPROPS, double *COORDS, int MCRD, int NNODE, double *U, double *DU, double *V, double *A,
         int JTYPE, double *TIME, double DTIME, int KSTEP, int KINC, int JELEM, double *PARAMS, int NDLOAD, int *JDLTYP,
         double *ADLMAG, double *PREDEF, int NPREDF, int *LFLAGS, int MLVARX, double *DDLMAG, int MDLOAD, double *PNEWDT,
         double *JPROPS, int NJPROP, double PERIOD) {

    // Define Gauss points
    double x_Gspoint[4] = {-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626};
    double y_Gspoint[4] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};

    // Initialize global variables
    if (KINC == 1) {
        ele_bond_state_pd[0][0][(int)COORDS[1 * NNODE + 8]][(int)COORDS[0 * NNODE + 8]] = 0;
        for (int i = 0; i < NSVARS; i++) {
            SVARS[i] = 0.0;
        }
        for (int i = 0; i < 4; i++) {
            ele_bond_sum_state_fem[i][0][(int)COORDS[0 * NNODE + 8]] = 0.0;
            ele_bond_sum_state_pd[i][0][(int)COORDS[1 * NNODE + 8]][(int)COORDS[0 * NNODE + 8]] = 0;
        }
    }

    // Initialize matrices
    double curr_gauss[2][8], ref_gauss[2][8], curr_coord[2][8], block_stiff[16][16], rhsmatrx[NDOFEL][1];
    mat_zero((double *)curr_gauss, 2, 8);
    mat_zero((double *)ref_gauss, 2, 8);
    mat_zero((double *)curr_coord, 2, 8);
    mat_zero((double *)AMATRX, NDOFEL, NDOFEL);
    mat_zero((double *)block_stiff, 16, 16);
    mat_zero((double *)rhsmatrx, NDOFEL, NRHS);

    // Calculate current node coordinates
    for (int i = 0; i < 8; i++) {
        curr_coord[0][i] = COORDS[0 * NNODE + i] + U[i * 2];
        curr_coord[1][i] = COORDS[1 * NNODE + i] + U[i * 2 + 1];
    }

    // Calculate Gauss coordinates
    gauss_loc(curr_coord, curr_gauss);
    gauss_loc((double (*)[8])COORDS, ref_gauss);

    // Determine bond breakage
    int loop_number = 1;
    int new_bre_bond = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double ref_bond = sqrt(pow(ref_gauss[0][i] - ref_gauss[0][j + 4], 2) + pow(ref_gauss[1][i] - ref_gauss[1][j + 4], 2));
            double curr_bond = sqrt(pow(curr_gauss[0][i] - curr_gauss[0][j + 4], 2) + pow(curr_gauss[1][i] - curr_gauss[1][j + 4], 2));
            if (SVARS[2 + loop_number] == 0.0) {
                if (ref_bond != 0.0) {
                    if (((curr_bond - ref_bond) / ref_bond) >= PROPS[1]) {
                        new_bre_bond++;
                        SVARS[2 + loop_number] = 1.0;
                    }
                }
            }
            loop_number++;
        }
    }

    // Assemble element stiffness matrix
    loop_number = 1;
    int para_build_bond = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (SVARS[2 + loop_number] != 1.0) {
                calc_gauss_intg(x_Gspoint[i], y_Gspoint[i], (double (*)[8])COORDS, x_Gspoint[j], y_Gspoint[j], block_stiff, PROPS[0], PROPS[2], PROPS[3]);
                for (int k = 0; k < NDOFEL; k++) {
                    for (int l = 0; l < NDOFEL; l++) {
                        AMATRX[k * NDOFEL + l] += block_stiff[k * 16 + l];
                    }
                }
                para_build_bond++;
            }
            loop_number++;
        }
    }

    mat_mult((double *)AMATRX, U, (double *)rhsmatrx, NDOFEL, NDOFEL, NDOFEL, 1);
    for (int i = 0; i < NDOFEL; i++) {
        RHS[i] = -rhsmatrx[i][0];
    }

    // Update state variables
    SVARS[0] = new_bre_bond;
    SVARS[19] = para_build_bond;
    SVARS[20] = 16 - para_build_bond;
    SVARS[1] = (16 - para_build_bond) / 16.0;

    // Store the number of broken bonds
    ele_bond_state_pd[0][0][(int)COORDS[1 * NNODE + 8]][(int)COORDS[0 * NNODE + 8]] = 1;
    for (int i = 0; i < 4; i++) {
        ele_bond_sum_state_pd[i][0][(int)COORDS[1 * NNODE + 8]][(int)COORDS[0 * NNODE + 8]] = SVARS[20];
    }
}
#endif


int main() {
    // Example usage of the functions
    double xRef[2][8] = {{0}};
    double loc_gauss[2][8] = {{0}};
    double block_stiff[16][16] = {{0}};
    double ele_esize = 1.0, param_1 = 1.0, param_2 = 1.0;

    gauss_loc(xRef, loc_gauss);
    calc_gauss_intg(0.0, 0.0, xRef, 0.0, 0.0, block_stiff, ele_esize, param_1, param_2);

    return 0;
}

