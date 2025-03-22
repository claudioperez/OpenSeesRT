//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <MatrixND.h>
#include <Matrix3D.h>
#include <Triad.h>
#include <Vector3D.h>
#include <Rotations.hpp>

class CrisfieldTransform {
public:
    CrisfieldTransform() {}

    int
    update(const Versor& qI, const Versor& qJ, const Vector3D& dx)
    {

        Ln = dx.norm();

        {
            Vector3D gammaw = CayleyFromVersor(qJ.mult_conj(qI));

            gammaw *= 0.5;

        //  Qbar = VersorProduct(VersorFromMatrix(CaySO3(gammaw)), qI);
            Qbar = VersorFromMatrix(CaySO3(gammaw)*MatrixFromVersor(qI));
            Triad r{CaySO3(gammaw)*MatrixFromVersor(qI)};
            r1 = r[1];
            r2 = r[2];
            r3 = r[3];
        }

        //
        // Compute the base vectors e2, e3
        //
        {
            // 'rotate' the mean rotation matrix Rbar on to e1 to
            // obtain e2 and e3 (using the 'mid-point' procedure)
            //
            // Vector3D e1, e2, e3;
            e[0]  = dx;
            e[0] /= Ln;
            Triad r = Triad{MatrixFromVersor(Qbar)};
            Vector3D r1 = r[1],
                     r2 = r[2],
                     r3 = r[3];

            // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);
        
            Vector3D tmp;
            tmp  = e[0];
            tmp += r1;//Qbar.rotate(E1);
        
            e[1] = tmp;
            {
                // const Vector3D r2 = Qbar.rotate(E2);
                e[1] *= 0.5*r2.dot(e[0]);
                e[1].addVector(-1.0,  r2, 1.0);
            }
        
            // e3 = r3 - (e1 + r1)*((r3^e1)*0.5);
            e[2] = tmp;
            {
                // const Vector3D r3 = Qbar.rotate(E3);
                e[2] *= r3.dot(e[0])*0.5;
                e[2].addVector(-1.0,  r3, 1.0);
            }
        }
        return 0;
    }

    constexpr const Vector3D& 
    getBasisE1() const noexcept
    {
      return e[0];
    }
    constexpr const Vector3D&
    getBasisE2() const noexcept
    {
      return e[1];
    }
    constexpr const Vector3D&
    getBasisE3() const noexcept
    {
      return e[2];
    }

    inline Matrix3D
    getRotation() const noexcept
    {
      Matrix3D E;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            E(i,j) = e[j][i];
      return E;
    }

    const Versor&
    getReference()
    {
      return Qbar;
    }
#if 0
    int
    addTangent(MatrixND<12,12>& kg, const VectorND<12>& pl)
    {    
        const Triad rI{MatrixFromVersor(Q_pres[0])},
                    rJ{MatrixFromVersor(Q_pres[1])};
        const Vector3D 
        &e1  =  getBasisE1(), // E[1],
        &e2  =  getBasisE2(), // E[2],
        &e3  =  getBasisE3(), // E[3],
        &rI1 = rI[1], // .rotate(E1), 
        &rI2 = rI[2], // .rotate(E2), 
        &rI3 = rI[3], // .rotate(E3),
        &rJ1 = rJ[1], // .rotate(E1), 
        &rJ2 = rJ[2], // .rotate(E2), 
        &rJ3 = rJ[3]; // .rotate(E3);
        // NOTE[cmp] 
        // CorotFrameTransf3d03::compTransfMatrixBasicGlobal must be 
        // called first to set Lr1, Lr2 and T

        // Matrix3D A;
        // for (int i = 0; i < 3; i++)
        //   for (int j = 0; j < 3; j++)
        //     A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;
        // getLMatrix(A, e1, r1, r2, Lr2);
        // getLMatrix(A, e1, r1, r3, Lr3);

        //
        // Ksigma1
        //
        {
        const double N = -pl[0]; // Axial force
        // a=0
        kg.assemble(A, 0, 0,  N);
        kg.assemble(A, 0, 6, -N);
        // a=1
        kg.assemble(A, 6, 0, -N);
        kg.assemble(A, 6, 6,  N);
        }

        //
        // Ksigma3
        //
        //  ks3 = [o kbar2  |  o kbar4];
        //
        //  where
        //
        //    kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;
        //
        //    kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));
        //
        // or
        //
        //  ks3 = [o ka+kb  |  o kc+kd];
        //      = [o ka     |  o kc] + [o kb  |  o kd];
        //
        //  where
        //
        //    ka = -Lr2*S(rI3)*m(3)  
        //         +Lr2*S(rI1)*m(1);
        //    kb =  Lr3*S(rI2)*m(3)  
        //         -Lr3*S(rI1)*m(2);
        //
        //    kc =  Lr2*S(rJ3)*m(3)
        //         -Lr2*S(rJ1)*m(4);
        //    kd = -Lr3*S(rJ2)*m(3)  
        //         +Lr3*S(rJ1)*m(5);

        static VectorND<6> m;
        m[0] =  0.5*pl[imx]/std::cos(ul(imx));
        m[2] = -0.5*pl[imy]/std::cos(ul(imy));
        m[1] =  0.5*pl[imz]/std::cos(ul(imz));

        m[3] =  0.5*pl[jmx]/std::cos(ul(jmx));
        m[5] = -0.5*pl[jmy]/std::cos(ul(jmy));
        m[4] =  0.5*pl[jmz]/std::cos(ul(jmz));


        static Matrix3D Sm;
        Sm.zero();
        Sm.addSpin(rI3,  m[3]);
        Sm.addSpin(rI1,  m[1]);
        static MatrixND<12,3> kbar;
        kbar.zero();
        kbar.addMatrixProduct(Lr2, Sm, -1.0);

        Sm.zero();
        Sm.addSpin(rI2,  m[3]);
        Sm.addSpin(rI1, -m[2]);
        kbar.addMatrixProduct(Lr3, Sm,  1.0);

        kg.assemble(kbar, 0, 3, 1.0);
        kg.assembleTranspose(kbar, 3, 0, 1.0);

        Sm.zero();
        Sm.addSpin(rJ3,  m[3]);
        Sm.addSpin(rJ1, -m[4]);
        kbar.zero();
        kbar.addMatrixProduct(Lr2, Sm, 1.0);

        Sm.zero();
        Sm.addSpin(rJ2, m[3]);
        Sm.addSpin(rJ1, m[5]);
        kbar.addMatrixProduct(Lr3, Sm,  -1.0);

        kg.assemble(kbar, 0, 9, 1.0);
        kg.assembleTranspose(kbar, 9, 0, 1.0);


        //
        // Ksigma4
        //
        {
        static Matrix3D ks33;
    
        ks33.zero();
        ks33.addSpinProduct(e2, rI3,  m[3]);
        ks33.addSpinProduct(e3, rI2, -m[3]);
        ks33.addSpinProduct(e2, rI1,  m[1]);
        ks33.addSpinProduct(e1, rI2, -m[1]);
        ks33.addSpinProduct(e3, rI1,  m[2]);
        ks33.addSpinProduct(e1, rI3, -m[2]);
        kg.assemble(ks33, 3, 3, 1.0);
        }

        //
        // Ksigma4
        //
        {
        static Matrix3D ks33;
        ks33.zero();
        ks33.addSpinProduct(e2, rJ3, -m[3]);
        ks33.addSpinProduct(e3, rJ2,  m[3]);
        ks33.addSpinProduct(e2, rJ1,  m[4]);
        ks33.addSpinProduct(e1, rJ2, -m[4]);
        ks33.addSpinProduct(e3, rJ1,  m[5]);
        ks33.addSpinProduct(e1, rJ3, -m[5]);
    
        kg.assemble(ks33, 9, 9, 1.0);
        }


        //
        // Ksigma5
        //
        //  Ks5 = [ Ks5_11   Ks5_12 | -Ks5_11   Ks5_14;
        //          Ks5_12'    O    | -Ks5_12'   O;
        //         -Ks5_11  -Ks5_12 |  Ks5_11  -Ks5_14;
        //          Ks5_14t     O   | -Ks5_14'   O];
        //
        //
        // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
        //   = 1/Ln * (m[1]*rI2 + m[2]*rI3)
        //   + 1/Ln * (m[4]*rJ2 + m[5]*rJ3);
        //   = vi + vj
        //
        {
        OPS_STATIC Vector3D v;
        v.addVector(0.0, rI2, m[1]);
        v.addVector(1.0, rI3, m[2]);
        v.addVector(1.0, rJ2, m[4]);
        v.addVector(1.0, rJ3, m[5]);
        v /= Ln;

        // Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
        //        = A*vi*e1' + e1*vi'*A + (e1'*vi)*A
        //        + A*vj*e1' + e1*vj'*A + (e1'*vj)*A;
        //
        Matrix3D ks33;
        ks33.zero();
        ks33.addMatrix(A, e1.dot(v));

        static Matrix3D m33;
        m33.zero();
        m33.addTensorProduct(v, e1, 1.0);

        ks33.addMatrixProduct(A, m33, 1.0);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            m33(i,j) = e1[i]*v[j];

        ks33.addMatrixProduct(m33, A, 1.0);

        kg.assemble(ks33, 0, 0,  1.0);
        kg.assemble(ks33, 0, 6, -1.0);
        kg.assemble(ks33, 6, 0, -1.0);
        kg.assemble(ks33, 6, 6,  1.0);
        }

        // Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));

        Matrix3D ks33;
        ks33.zero();
        ks33.addMatrixSpinProduct(A, rI2, -m[1]);
        ks33.addMatrixSpinProduct(A, rI3, -m[2]);

        kg.assemble(ks33, 0, 3,  1.0);
        kg.assemble(ks33, 6, 3, -1.0);
        kg.assembleTranspose(ks33, 3, 0,  1.0);
        kg.assembleTranspose(ks33, 3, 6, -1.0);

        //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));

        ks33.zero();
        ks33.addMatrixSpinProduct(A, rJ2, -m[4]);
        ks33.addMatrixSpinProduct(A, rJ3, -m[5]);

        kg.assemble(ks33, 0, 9,  1.0);
        kg.assemble(ks33, 6, 9, -1.0);

        kg.assembleTranspose(ks33, 9, 0,  1.0);
        kg.assembleTranspose(ks33, 9, 6, -1.0);

        // Ksigma -------------------------------
        OPS_STATIC Vector3D rm;

        rm = rI3;
        rm.addVector(1.0, rJ3, -1.0);
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rm), m[3]);

    //  rm = rJ2;
        rm.addVector(0.0, rJ2, -1.0);
        rm.addVector(1.0, rI2, -1.0);
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3,  rm), m[3]);
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rI1), m[1]);
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3, rI1), m[2]);
        //
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rJ1), m[4]);
        kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3, rJ1), m[5]);

        //
        //  T' * diag (M .* tan(thetal))*T
        //

        for (int node=0; node<2; node++) {
          for (int k = 0; k < 3; k++) {
              const double factor =  pl[6*node+3+k] * std::tan(ul[(node ? jmx : imx) + k]);
              for (int i = 0; i < 12; i++) {
              const double Tki = T((node ? jmx : imx) + k,i);
              for (int j = 0; j < 12; j++)
                  kg(i,j) += Tki * factor * T((node ? jmx : imx) + k, j);
              }
          }
        }

        return 0;
    }
#endif
private:
    Versor Qbar;
    Vector3D r1, r2, r3;
    Vector3D e[3];
    double Ln;

    // Auxiliary
//  Matrix3D A;
//  MatrixND<12,3> Lr2, Lr3;   // auxiliary matrices
};
