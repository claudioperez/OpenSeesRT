//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Implementation of a 3D Vector
//
// Original implementation: Massimo Petracca (ASDEA)
//
#ifndef Vector3D_h
#define Vector3D_h

#include <VectorND.h>

#ifndef M_PI
#  define M_PI 3.1415926535897932384626433832795
#endif


namespace OpenSees {class Matrix3D;};

using Vector3D = OpenSees::VectorND<3, double>;

#if 0
// template<class scalar_t=double>
class Vector3D: public OpenSees::VectorND<3>
{
public:
    typedef double scalar_t;

    friend class OpenSees::Matrix3D;

    /**
     * Creates a Zero Vector3D.
    */
//  Vector3D()
//  {
//      values[0] = values[1] = values[2] = 0.0;
//  }

    /**
    Creates a Vector3D from its coefficients.
    @param x x coefficient
    @param y y coefficient
    @param z z coefficient
    */
//  Vector3D(scalar_t x, scalar_t y, scalar_t z)
//  {
//      values[0] = x;
//      values[1] = y;
//      values[2] = z;
//  }

    /**
    Creates a Vector3D from a Vector starting from index.
    Note: no check is made on the input data size!
    @param v the vector
    @param index the index (optional, default = 0)
    */
//  Vector3D(const Vector& v, size_t index = 0)
//  {
//      values[0] = v(0 + index);
//      values[1] = v(1 + index);
//      values[2] = v(2 + index);
//  }

    /**
    Creates a Vector3D from another Vector3D.
    @param other the other Vector3D
    */
//  Vector3D(const Vector3D& other)
//  {
//      values[0] = other.values[0];
//      values[1] = other.values[1];
//      values[2] = other.values[2];
//  }

public:

//  operator Vector() { 
//    return Vector(values, 3);
//  }

    /**
    Copies a Vector3D.
    @param other the other Vector3D
    */
    inline Vector3D& 
    operator= (const Vector3D& other)
    {
        if (this != &other) {
            values[0] = other.values[0];
            values[1] = other.values[1];
            values[2] = other.values[2];
        }
        return *this;
    }

    template<typename VecT> inline
    Vector3D &operator=(const VecT &right) {
      for (int i=0; i< 3; i++)
        values[i] = right[i];
      return *this;
    }

    template <class VecT> inline
    Vector3D &operator-=(const VecT &right) {
      for (int i=0; i< 3; i++)
        values[i] -= right[i];
      return *this;
    }

    template <class VecT> inline
    Vector3D &operator+=(const VecT &right) {
      for (int i=0; i< 3; i++)
        values[i] += right[i];
      return *this;
    }

public:

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline constexpr scalar_t operator()(size_t i) const { return values[i]; }

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline constexpr scalar_t& operator()(size_t i) { return values[i]; }

public:

    /**
    Returns the squared norm this vector.
    @return the squared norm of this vector.
    */
    inline scalar_t squaredNorm() const
    {
        return 
            values[0] * values[0] +
            values[1] * values[1] +
            values[2] * values[2];
    }

    /**
    makes this vector a unit vector.
    @return the norm of this vector.
    */
    inline scalar_t normalize() {
        scalar_t n = norm();
        if (n > 0.0) {
            values[0] /= n;
            values[1] /= n;
            values[2] /= n;
        }
        return n;
    }

    /**
    Returns the dot product this vector with another vector.
    @param b the other vector
    @return the dot product.
    */
//  inline constexpr scalar_t 
//  dot(const Vector3D& b) const {
//      return
//          values[0] * b.values[0] +
//          values[1] * b.values[1] +
//          values[2] * b.values[2];
//  }


    /**
    Returns the cross product this vector with another vector.
    @param b the other vector
    @return the cross product.
    */
//  template <class VecB, class VecC> inline 
//  void cross(const VecB& b, VecC& c) const {
//      c[0] = values[1] * b[2] - values[2] * b[1];
//      c[1] = values[2] * b[0] - values[0] * b[2];
//      c[2] = values[0] * b[1] - values[1] * b[0];
//      return c;
//  }


//  template <class Vec3T> inline
//  Vector3D cross(const Vec3T& b) const {
//      Vector3D c;
//      c.values[0] = values[1] * b.values[2] - values[2] * b.values[1];
//      c.values[1] = values[2] * b.values[0] - values[0] * b.values[2];
//      c.values[2] = values[0] * b.values[1] - values[1] * b.values[0];
//      return c;
//  }

#if 0
    void addCrossProduct(const Vector3D& a, const Vector3D& b, const double scale) {
        values[0] += scale*(a.values[1] * b.values[2] - a.values[2] * b.values[1]);
        values[1] += scale*(a.values[2] * b.values[0] - a.values[0] * b.values[2]);
        values[2] += scale*(a.values[0] * b.values[1] - a.values[1] * b.values[0]);
    }
#endif

public:

//  inline void operator += (const Vector3D& b) {
//      values[0] += b.values[0];
//      values[1] += b.values[1];
//      values[2] += b.values[2];
//  }

//  inline Vector3D operator + (const Vector3D& b) const {
//      Vector3D a(*this);
//      a += b;
//      return a;
//  }

//  inline void operator -= (const Vector3D& b) {
//      values[0] -= b.values[0];
//      values[1] -= b.values[1];
//      values[2] -= b.values[2];
//  }

//  inline Vector3D operator - (const Vector3D& b) const {
//      Vector3D a(*this);
//      a -= b;
//      return a;
//  }

//  inline void operator *= (scalar_t b) {
//      values[0] *= b;
//      values[1] *= b;
//      values[2] *= b;
//  }

//  inline Vector3D operator * (scalar_t b) const {
//      Vector3D a(*this);
//      a *= b;
//      return a;
//  }

//  inline void operator /= (scalar_t b) {
//      values[0] /= b;
//      values[1] /= b;
//      values[2] /= b;
//  }

//  inline Vector3D operator / (scalar_t b) const {
//      Vector3D a(*this);
//      a /= b;
//      return a;
//  }
};


inline Vector3D operator * (double a, const Vector3D& b) {
    return b * a;
}
#endif


#endif // Vector3D_h
