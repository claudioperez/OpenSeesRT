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

    const Matrix3D&
    getRotation()
    {
      return e;
    }

    const Versor&
    getReference()
    {
      return Qbar;
    }

    int
    update(const Versor& qI, const Versor& qJ, const Vector3D& dx)
    {

        Ln = dx.norm();

        {
            Vector3D gammaw = CayleyFromVersor(qJ.mult_conj(qI));

            gammaw *= 0.5;

        //  Qbar = VersorProduct(VersorFromMatrix(CaySO3(gammaw)), qI);
            Qbar = VersorFromMatrix(CaySO3(gammaw)*MatrixFromVersor(qI));
        }

        //
        // Compute the base vectors e2, e3
        //
        {
            // 'rotate' the mean rotation matrix Rbar on to e1 to
            // obtain e2 and e3 (using the 'mid-point' procedure)
            //
            Vector3D e1, e2, e3;
            e1  = dx;
            e1 /= Ln;
            Triad r = Triad{MatrixFromVersor(Qbar)};
            Vector3D r1 = r[1],
                     r2 = r[2],
                     r3 = r[3];

            // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);
        
            Vector3D tmp;
            tmp  = e1;
            tmp += r1;//Qbar.rotate(E1);
        
            e2 = tmp;
            {
                // const Vector3D r2 = Qbar.rotate(E2);
                e2 *= 0.5*r2.dot(e1);
                e2.addVector(-1.0,  r2, 1.0);
            }
        
            // e3 = r3 - (e1 + r1)*((r3^e1)*0.5);
            e3 = tmp;
            {
                // const Vector3D r3 = Qbar.rotate(E3);
                e3 *= r3.dot(e1)*0.5;
                e3.addVector(-1.0,  r3, 1.0);
            }
            for (int k = 0; k < 3; k ++) {
                e(k,0) = e1[k];
                e(k,1) = e2[k];
                e(k,2) = e3[k];
            }
        }
        return 0;
    }

private:
    Versor Qbar;
    Matrix3D e;
    double Ln;
    constexpr static Vector3D E1 {1, 0, 0}, 
                              E2 {0, 1, 0},
                              E3 {0, 0, 1};
};