#pragma once
#include <cmath>
#include <VectorND.h>

struct Versor: public OpenSees::VectorND<4> {
//  inline VectorND<3>
//  as_vector()
//  {
//    VectorND<3> v;
//  }

    template<typename Vec3T>
    static inline Versor
    from_vector(const Vec3T  &theta)
    {

        double t = 0.0; // theta.norm();
        for (int i=0; i<3; i++)
            t += theta[i]*theta[i];
        t = std::sqrt(t);

        Versor q;
        if (t == 0)
            q.zero();

        else {
            const double factor = std::sin(t*0.5) / t;
            for (int i = 0; i < 3; i++)
                q[i] = theta[i] * factor;
        }

        // Scalar part
        q[3] = std::cos(t*0.5);

        return q;
    }

  inline Versor &operator*= (const Versor& b)
  {

      values[0] = values[3]*b[0] + values[0]*b[3] + values[1]*b[2] - values[2]*b[1];
      values[1] = values[3]*b[1] + values[1]*b[3] + values[2]*b[0] - values[0]*b[2];
      values[2] = values[3]*b[2] + values[2]*b[3] + values[0]*b[1] - values[1]*b[0];
      values[3] = values[3]*b[3] - values[0]*b[0] - values[1]*b[1] - values[2]*b[2];
      return *this;
  }
};

#if 1
inline Versor
operator*(const Versor &qa, const Versor &qb)
{
    const double qa0 = qa[0],
                 qa1 = qa[1],
                 qa2 = qa[2],
                 qa3 = qa[3],
                 qb0 = qb[0],
                 qb1 = qb[1],
                 qb2 = qb[2],
                 qb3 = qb[3];

    // Calculate the dot product qa.qb
    const double qaTqb = qa0*qb0 + qa1*qb1 + qa2*qb2;

    // Calculate the cross-product qa x qb
    const double
      qaxqb0 = qa1*qb2 - qa2*qb1,
      qaxqb1 = qa2*qb0 - qa0*qb2,
      qaxqb2 = qa0*qb1 - qa1*qb0;

    // Calculate the quaternion product
    Versor q12;
    q12[0] = qa3*qb0 + qb3*qa0 - qaxqb0;
    q12[1] = qa3*qb1 + qb3*qa1 - qaxqb1;
    q12[2] = qa3*qb2 + qb3*qa2 - qaxqb2;
    q12[3] = qa3*qb3 - qaTqb;
    return q12;
}
#else
inline Versor 
operator* (const Versor& a, const Versor& b)
{
    return Versor{{
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] + a[2]*b[0] + a[3]*b[1] - a[1]*b[3],
        a[0]*b[3] + a[3]*b[0] + a[1]*b[2] - a[2]*b[1]
    }};
}
#endif