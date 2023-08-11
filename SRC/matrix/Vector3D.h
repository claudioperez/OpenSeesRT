/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Implementation of a 3D Vector
//
// Original implementation: Massimo Petracca (ASDEA)
//
#ifndef Vector3D_h
#define Vector3D_h

#include <math.h>
#include <stdint.h>
#include <cmath>
#include <limits>
#include <VectorND.h>
#include <Vector.h>
#include <Matrix.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI


template<class scalar_t=double>
class Vector3D: public OpenSees::VectorND<3>
{
public:

    /**
     * Creates a Zero Vector3D.
    */
    Vector3D()
    {
        mData[0] = mData[1] = mData[2] = 0.0;
    }

    /**
    Creates a Vector3D from its coefficients.
    @param x x coefficient
    @param y y coefficient
    @param z z coefficient
    */
    Vector3D(scalar_t x, scalar_t y, scalar_t z)
    {
        mData[0] = x;
        mData[1] = y;
        mData[2] = z;
    }

    /**
    Creates a Vector3D from a Vector starting from index.
    Note: no check is made on the input data size!
    @param v the vector
    @param index the index (optional, default = 0)
    */
    Vector3D(const Vector& v, size_t index = 0)
    {
        mData[0] = v(0 + index);
        mData[1] = v(1 + index);
        mData[2] = v(2 + index);
    }

    /**
    Creates a Vector3D from another Vector3D.
    @param other the other Vector3D
    */
    Vector3D(const Vector3D& other)
    {
        mData[0] = other.mData[0];
        mData[1] = other.mData[1];
        mData[2] = other.mData[2];
    }

public:

    /**
    Copies a Vector3D.
    @param other the other Vector3D
    */
    Vector3D& operator= (const Vector3D& other)
    {
        if (this != &other) {
            mData[0] = other.mData[0];
            mData[1] = other.mData[1];
            mData[2] = other.mData[2];
        }
        return *this;
    }

public:

    /**
    Returns the X coefficient of this vector.
    @return the X coefficient of this vector.
    */
    inline const scalar_t x()const { return mData[0]; }

    /**
    Returns the Y coefficient of this vector.
    @return the Y coefficient of this vector.
    */
    inline const scalar_t y()const { return mData[1]; }

    /**
    Returns the Z coefficient of this vector.
    @return the Z coefficient of this vector.
    */
    inline const scalar_t z()const { return mData[2]; }

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline constexpr scalar_t operator()(size_t i) const { return mData[i]; }

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline constexpr scalar_t& operator()(size_t i) { return mData[i]; }

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline scalar_t operator[](size_t i) const { return mData[i]; }

    /**
    Returns the i-th coefficient of this vector.
    @return the i-th coefficient of this vector.
    */
    inline scalar_t& operator[](size_t i) { return mData[i]; }

public:

    /**
    Returns the squared norm this vector.
    @return the squared norm of this vector.
    */
    inline scalar_t squaredNorm() const
    {
        return 
            mData[0] * mData[0] +
            mData[1] * mData[1] +
            mData[2] * mData[2];
    }

    /**
    Returns the norm this vector.
    @return the norm of this vector.
    */
    inline scalar_t norm() const
    {
        return std::sqrt(squaredNorm());
    }

    /**
    makes this vector a unit vector.
    @return the norm of this vector.
    */
    inline scalar_t normalize() {
        scalar_t n = norm();
        if (n > 0.0) {
            mData[0] /= n;
            mData[1] /= n;
            mData[2] /= n;
        }
        return n;
    }

    /**
    Returns the dot product this vector with another vector.
    @param b the other vector
    @return the dot product.
    */
    inline constexpr scalar_t dot(const Vector3D& b) const {
        return
            mData[0] * b.mData[0] +
            mData[1] * b.mData[1] +
            mData[2] * b.mData[2];
    }

    /**
    Returns the cross product this vector with another vector.
    @param b the other vector
    @return the cross product.
    */
    inline constexpr Vector3D cross(const Vector3D& b) const {
        Vector3D c;
        c.mData[0] = mData[1] * b.mData[2] - mData[2] * b.mData[1];
        c.mData[1] = mData[2] * b.mData[0] - mData[0] * b.mData[2];
        c.mData[2] = mData[0] * b.mData[1] - mData[1] * b.mData[0];
        return c;
    }

public:

    inline void operator += (const Vector3D& b) {
        mData[0] += b.mData[0];
        mData[1] += b.mData[1];
        mData[2] += b.mData[2];
    }

    inline Vector3D operator + (const Vector3D& b) const {
        Vector3D a(*this);
        a += b;
        return a;
    }

    inline void operator -= (const Vector3D& b) {
        mData[0] -= b.mData[0];
        mData[1] -= b.mData[1];
        mData[2] -= b.mData[2];
    }

    inline Vector3D operator - (const Vector3D& b) const {
        Vector3D a(*this);
        a -= b;
        return a;
    }

    inline void operator *= (scalar_t b) {
        mData[0] *= b;
        mData[1] *= b;
        mData[2] *= b;
    }

    inline Vector3D operator * (scalar_t b) const {
        Vector3D a(*this);
        a *= b;
        return a;
    }

    inline void operator /= (scalar_t b) {
        mData[0] /= b;
        mData[1] /= b;
        mData[2] /= b;
    }

    inline Vector3D operator / (scalar_t b) const {
        Vector3D a(*this);
        a /= b;
        return a;
    }

private:
    scalar_t mData[3];
};

template<class scalar_t>
inline Vector3D<scalar_t> operator * (scalar_t a, const Vector3D<scalar_t>& b) {
    return b * a;
}

/**
Prints this vector to a input stream
@param s the output stream
@param v the vector
@return the stream
*/
template<class TStream, class T>
inline TStream& operator << (TStream& s, const Vector3D<T>& v) {
    return (s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")");
}

#endif // Vector3D_h
