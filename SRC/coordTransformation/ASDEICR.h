/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of the EICR (Element Independent CoRotational) formulation
// 
// References:
//
// Felippa, Carlos A. 
//   "A systematic approach to the element-independent corotational 
//   dynamics of finite elements". 
//   Technical Report CU-CAS-00-03, Center for Aerospace Structures, 2000.
//
// Felippa, Carlos A., and Bjorn Haugen. 
//   "A unified formulation of small-strain corotational finite elements: I. Theory." 
//   Computer Methods in Applied Mechanics and Engineering 194.21-24 (2005): 2285-2335.
//
#ifndef ASDEICR_h
#define ASDEICR_h

#include <Quaternion.h>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <vector>

/** EICR Element Independent Corotational formulation
*
* E.I.C.R. is a utility class containing static methods related to
* the Element Independent Corotational Formulation.
*/
namespace EICR
{

typedef std::vector<Vector3D> NodeContainerType;

typedef ASDQuaternion<double> QuaternionType;

using OpenSees::Matrix3D;
using std::size_t;

/**
* Set the input matrix to the zero matrix of requested size. Resize is done if necessary
* @param n the number of rows
* @param n the number of columns
* @param I the input/output matrix
*/
inline static void SetZero(size_t n, size_t m, Matrix& I)
{
    if ((I.noRows() != n) || (I.noCols() != m))
        I.resize(n, m);
    I.Zero();
}

/**
* Sets the input matrix to the identity of requested size. Resize is done if necessary
* @param n the number of rows and columns
* @param I the input/output matrix
* @param value (optional, default = 1) the value on the diagonal terms
*/
inline static void SetIdentity(size_t n, Matrix& I, double value = 1.0)
{
    SetZero(n, n, I);
    for (size_t j = 0; j < n; ++j)
        I(j, j) = value;
}

/**
* Copy the block [begin:end[ from A to B. B size must be end-begin. 
*
* @param A the first vector
* @param begin the first index
* @param end the last+1 index
* @param B the second vector
*/
template <typename VectType>
inline static void GetBlock(const Vector& A, size_t begin, size_t end, VectType& B)
{
    size_t n = end - begin;
    for (size_t i = 0; i < n; ++i)
        B(i) = A(i + begin);
}

/**
* Copy the block [begin:end[ from A to B. B size must be end-begin.
*
* @param A the first matrix
* @param begin the first index
* @param end the last+1 index
* @param B the second matrix
*/
template <typename MatT>
inline static void GetBlock(const Matrix& A, size_t begin, size_t end, MatT& B)
{
    size_t n = end - begin;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            B(i, j) = A(i + begin, j + begin);
}

/**
* Copy the block [begin:end[ from B to A. B size must be end-begin.
*
* @param A the first vector
* @param begin the first index
* @param end the last+1 index
* @param B the second vector
*/
template<class MatT>
inline static void SetBlock(Matrix& A, size_t begin, size_t end, const MatT& B)
{
    size_t n = end - begin;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            A(i + begin, j + begin) = B(i, j);
}

/**
* Computes the Translational Projector Matrix.
* The output is a square matrix of size num_nodes*6.
* Note that 6 Degrees Of Freedom are assumed for each node.
* @param num_nodes the number of nodes
* @param P the Translational Projector Matrix
*/
inline static void Compute_Pt(size_t num_nodes, Matrix& P)
{
    double a = double(num_nodes - 1) / double(num_nodes);
    double b = -1.0 / double(num_nodes);

    size_t num_dofs = num_nodes * 6;

    SetIdentity(num_dofs, P);

    for (size_t i = 0; i < num_nodes; i++) {
        size_t j = i * 6;

        // diagonal block
        P(j, j) = a;
        P(j + 1, j + 1) = a;
        P(j + 2, j + 2) = a;

        // out-of-diagonal block
        for (size_t k = i + 1; k < num_nodes; k++) {
            size_t w = k * 6;

            P(j, w) = b;
            P(j + 1, w + 1) = b;
            P(j + 2, w + 2) = b;

            P(w, j) = b;
            P(w + 1, j + 1) = b;
            P(w + 2, j + 2) = b;
        }
    }
}

/**
* Computes the Spin Lever Matrix.
* The output is a rectangular matrix of 3 columns and nodes.size()*6 rows.
* Note that 6 Degrees Of Freedom are assumed for each node.
* @param nodes the input nodes
* @param S the Spin Lever Matrix
*/
inline static void Compute_S(const NodeContainerType& nodes, Matrix& S)
{
    size_t num_nodes = nodes.size();
    size_t num_dofs = num_nodes * 6;

    SetZero(num_dofs, 3, S);

    for (size_t i = 0; i < num_nodes; i++)
    {
        size_t j = i * 6;

        S.addSpinAtRow(nodes[i], -1.0, 0, j);

        S(j + 3, 0) = 1.0;
        S(j + 4, 1) = 1.0;
        S(j + 5, 2) = 1.0;
    }
}

/**
* Computes the Axial Vector Jacobian.
* The output is a square matrix of size displacements.size() (which is num_nodes * 6).
* Note that 6 Degrees Of Freedom are assumed for each node.
* @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
* @return the H matrix
*/
inline static void Compute_H(const Vector& displacements, Matrix& H)
{
    size_t num_dofs = displacements.Size();
    size_t num_nodes = num_dofs / 6;

    SetIdentity(num_dofs, H);

    static Matrix3D Hi;
    static Vector3D rv;

    for (size_t i = 0; i < num_nodes; i++) {
        size_t index = i * 6;

        GetBlock(displacements, index + 3, index + 6, rv);

        double angle = rv.norm();

        if (angle >= 2.0 * M_PI)
            angle = std::fmod(angle, 2.0 * M_PI);

        double eta;
        if (angle < 0.05) {
            double angle2 = angle  * angle;
            double angle4 = angle2 * angle2;
            double angle6 = angle4 * angle2;
            eta = 1.0 / 12.0 + 1.0 / 720.0 * angle2 + 1.0 / 30240.0 * angle4 + 1.0 / 1209600.0 * angle6;
        }
        else {
            eta = (1.0 - 0.5 * angle * std::tan(0.5 * M_PI - 0.5 * angle)) / (angle * angle);
        }

        // SO3::Spin(rv, Omega);
        // Omega2.addMatrixProduct(0.0, Omega, Omega, 1.0);

        // Hi = I - 0.5*Omega + eta*Omega*Omega
        Hi.zero();
        Hi.addDiagonal(1.0);
        Hi.addSpin(rv, -0.5);
        Hi.addSpinSquare(rv, eta);

        SetBlock(H, index + 3, index + 6, Hi);
    }
}

/**
* Compute the Spin derivative of (Axial Vector Jacobian)^T contracted with the nodal moment vector.
* The output is a square matrix of size displacements.size() (which is num_nodes * 6).
* Note that 6 Degrees Of Freedom are assumed for each node.
* @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
* @param forces the vector of nodal forces and moments in the local corotational coordinate system. (assumed size = num_nodes*6)
* @param H the Axial Vector Jacobian Matrix computed with a previous call to EICR::Compute_H(displacements)
* @return the L matrix
*/
inline static void Compute_L(const Vector& displacements, const Vector& forces, const Matrix& H, Matrix& L)
{
    size_t num_dofs = displacements.Size();
    size_t num_nodes = num_dofs / 6;

    SetZero(num_dofs, num_dofs, L);


    for (size_t i = 0; i < num_nodes; i++) {
        size_t index = i * 6;

        Vector3D rotationVector{ 0.0, 0.0, 0.0 };
        Vector3D momentVector{ 0.0, 0.0, 0.0 };
        GetBlock(displacements, index + 3, index + 6, rotationVector);
        GetBlock(forces,        index + 3, index + 6, momentVector);

        double angle = rotationVector.norm();

        if (angle >= 2.0 * M_PI)
            angle = std::fmod(angle, 2.0 * M_PI);

        double angle2 = angle  * angle;
        double angle4 = angle2 * angle2;
        double angle6 = angle4 * angle2;

        double mu;
        double eta;
        if (angle < 0.05) {
            eta = 1.0 / 12.0 + angle2 / 720.0 + angle4 / 30240.0 + angle6 / 1209600.0;
            mu = 1.0 / 360.0 + angle2 / 7560.0 + angle4 / 201600.0 + angle6 / 5987520.0;
        }
        else {
            eta = (1.0 - 0.5 * angle * std::tan(0.5 * M_PI - 0.5 * angle)) / (angle * angle);
            double sin_h_angle = std::sin(0.5 * angle);
            mu = (angle2 + 4.0 * std::cos(angle) + angle * std::sin(angle) - 4.0) / (4.0 * angle4 * sin_h_angle * sin_h_angle);
        }

        Matrix3D A{};
        A.addDiagonal(   rotationVector.dot(momentVector)*eta);
        A.addTensorProduct(rotationVector, momentVector,  eta);
        A.addTensorProduct(momentVector, rotationVector, -eta);

        Matrix3D Omega2{};
        Omega2.addSpinSquare(rotationVector, 1.0);

        A.addTensorProduct(Omega2*momentVector, rotationVector, mu);
        A.addSpin(momentVector, -0.5);

        Matrix3D Hi{};
        GetBlock(H, index + 3, index + 6, Hi);

        Matrix3D Li = A*Hi;
        SetBlock(L, index + 3, index + 6, Li);
    }
}
} // namespace EICR

#endif // !ASDEICR_h
