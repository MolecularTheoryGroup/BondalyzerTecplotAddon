#pragma once

#include "TASSERT.h"
#include "SimpleTypes.h"
#include <algorithm>

namespace tpcsm {

class Vec3 // for XYZ or UVW among other things
{
private:
    double m_x;
    double m_y;
    double m_z;
public:
    Vec3()
    {
        // leave uninitialized
    }
    Vec3(Vec3 const& other)
    {
        m_x = other.m_x;
        m_y = other.m_y;
        m_z = other.m_z;
    }
    Vec3(double xx, double yy, double zz)
    {
        m_x = xx;
        m_y = yy;
        m_z = zz;
    }
    double x() const { return m_x; }
    double y() const { return m_y; }
    double z() const { return m_z; }
    double& x() { return m_x; }
    double& y() { return m_y; }
    double& z() { return m_z; }

	// comparison

	bool operator ==(Vec3 const& other) const
	{
		return (m_x == other.m_x && m_y == other.m_y && m_z == other.m_z);
	}

    // unary minus
    Vec3 operator -() const
    {
        return Vec3(-m_x, -m_y, -m_z);
    }

    // note for vec3 not greater-than does not imply less-than-or-equal-to
    bool operator >(double val) const
    {
        return (m_x > val && m_y > val && m_z > val);
    }
    bool operator >=(double val) const
    {
        return (m_x >= val && m_y >= val && m_z >= val);
    }
    bool operator <(double val) const
    {
        return (m_x < val && m_y < val && m_z < val);
    }
    bool operator <=(double val) const
    {
        return (m_x <= val && m_y <= val && m_z < val);
    }

    // scalar addition
    Vec3& operator +=(double scalar)
    {
        m_x += scalar;
        m_y += scalar;
        m_z += scalar;
        return *this;
    }
    Vec3 operator +(double scalar) const
    {
        Vec3 result(*this);
        result += scalar;
        return result;
    }
    // member-wise addition
    Vec3& operator +=(Vec3 const& other)
    {
        m_x += other.m_x;
        m_y += other.m_y;
        m_z += other.m_z;
        return *this;
    }
    Vec3 operator +(Vec3 const& other) const
    {
        Vec3 result(*this);
        result += other;
        return result;
    }

    // scalar multiplication
    Vec3& operator *=(double scalar)
    {
        m_x *= scalar;
        m_y *= scalar;
        m_z *= scalar;
        return *this;
    }
    Vec3 operator *(double scalar) const
    {
        Vec3 result(*this);
        result *= scalar;
        return result;
    }
    // member-wise multiplication
    Vec3& operator *=(Vec3 const& other)
    {
        m_x *= other.m_x;
        m_y *= other.m_y;
        m_z *= other.m_z;
        return *this;
    }
    Vec3 operator *(Vec3 const& other) const
    {
        Vec3 result(*this);
        result *= other;
        return result;
    }

    // scalar subtraction
    Vec3& operator -=(double scalar)
    {
        *this = *this - scalar;
        return *this;
    }
    Vec3 operator -(double scalar) const
    {
        Vec3 result(
            m_x - scalar,
            m_y - scalar,
            m_z - scalar);
        return result;
    }
    // member-wise subtraction
    Vec3& operator -=(Vec3 const& other)
    {
        *this = *this - other;
        return *this;
    }
    Vec3 operator -(Vec3 const& other) const
    {
        Vec3 result(
            m_x - other.m_x,
            m_y - other.m_y,
            m_z - other.m_z);
        return result;
    }

    // scalar division
    Vec3& operator /=(double scalar)
    {
        *this = *this / scalar;
        return *this;
    }
    Vec3 operator /(double scalar) const
    {
        Vec3 result(
            m_x / scalar,
            m_y / scalar,
            m_z / scalar);
        return result;
    }
    // member-wise subtraction
    Vec3& operator /=(Vec3 const& other)
    {
        m_x /= other.m_x;
        m_y /= other.m_y;
        m_z /= other.m_z;
        return *this;
    }
    Vec3 operator /(Vec3 const& other) const
    {
        Vec3 result(*this);
        result /= other;
        return result;
    }

    // dot product
    inline double dot(Vec3 const& other) const
    {
        double dotProduct = m_x * other.m_x +
            m_y * other.m_y +
            m_z * other.m_z;
        return dotProduct;
    }

    // cross product
    inline Vec3 cross(Vec3 const& other) const
    {
        double const u1 = m_x;
        double const u2 = m_y;
        double const u3 = m_z;
        double const v1 = other.m_x;
        double const v2 = other.m_y;
        double const v3 = other.m_z;

        return Vec3(u2*v3 - u3*v2,
            u3*v1 - u1*v3,
            u1*v2 - u2*v1);
    }

    inline double getNorm() const
    {
        double const norm = sqrt(getNormSquared());
        return norm;
    }
    // save the sqrt call if not needed
    inline double getNormSquared() const
    {
        double const normSquared = this->dot(*this);
        return normSquared;
    }

    inline Vec3 normalize() const
    {
        double const norm = getNorm();
        if (norm > 0.0)
            return Vec3(*this / norm);
        else
            return *this;
    }

    static inline double epsilon() // TODO use constexpr when available (not in VS2013)
    {
        return 1.0e-02; // because vectors coming in from files might be floats we have to use a float epsilon
    }

    inline bool isNormalized() const
    {
        double const ns = getNormSquared();
        return std::abs(ns - 1.0) < epsilon();
    }

#ifndef NO_ASSERTS
    //
    inline bool nearlyZero() const
    {
        double const ns = getNormSquared();
        return ns < epsilon();
    }

    //
    inline bool nearlyEquals(Vec3 const& other) const
    {
        double dp = this->dot(other);
        return std::abs(dp - 1.0) < epsilon();
    }
#endif

    //
    inline Vec3 min(Vec3 const& other) const
    {
        return Vec3(std::min(m_x, other.m_x), std::min(m_y, other.m_y), std::min(m_z, other.m_z));
    }

    //
    inline Vec3 max(Vec3 const& other) const
    {
        return Vec3(std::max(m_x, other.m_x), std::max(m_y, other.m_y), std::max(m_z, other.m_z));
    }

    //
    static inline double tripleProduct(
        Vec3 const& a,
        Vec3 const& b,
        Vec3 const& c)
    {
        return a.cross(b).dot(c);
    }
};

// right-handed scalar addition
inline Vec3 operator +(double scalar, Vec3 const& vec3)
{
    Vec3 result(vec3);
    result += scalar;
    return result;
}
// right-handed scalar multiply
inline Vec3 operator *(double scalar, Vec3 const& vec3)
{
    Vec3 result(vec3);
    result *= scalar;
    return result;
}

}
