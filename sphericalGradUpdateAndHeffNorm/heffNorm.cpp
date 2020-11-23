#include "heffNorm.h"
#include "vec3.h"

// heffNorm.cpp
namespace tpcsm {
    //
    // Local class for Quaternions with unit norm are rotations in 3D space
    //
    namespace {
        class UnitQuaternion
        {
        private:
            double m_r; // There are two ways to represent Quaternions (r,i,j,k) and (r,V) where V is a 3D vector
            Vec3   m_v; // The (r,V) version makes the math simpler.

        private:
            UnitQuaternion( // this constructor is private mainly to prevent confusion with the (axis,angle) one
                double r,
                Vec3 const& v)
                : m_r(r)
                , m_v(v)
            {
                REQUIRE(isValid()); // only allowed to create unit quaternions
            }
            // is this constructor necessary?
            UnitQuaternion(
                double rVal,
                double iVal,
                double jVal,
                double kVal)
                : m_r(rVal)
                , m_v(iVal, jVal, kVal)
            {
                REQUIRE(isValid()); // only allowed to create unit quaternions
            }

        public:
            UnitQuaternion() // the identity rotation
                : m_r(1.0)
                , m_v(0.0, 0.0, 0.0)
            {}
            UnitQuaternion( // create quaternion that is a rotation of "angle" around vector "axis"
                Vec3 const& axis,
                double angle)
            {
                REQUIRE(axis.isNormalized());

                double const halfAngle = 0.5*angle;
                m_r = cos(halfAngle);
                m_v = axis * sin(halfAngle);

                ENSURE(isValid());
            }
        private:
            inline double normSquared() const // for UnitQuaternions this should always return 1, so gettting the norm is private to prevent abuse
            {
                double const ns = m_r*m_r + m_v.x()*m_v.x() + m_v.y()*m_v.y() + m_v.z()*m_v.z(); // norm squared
                return ns;
            }
            inline bool isNormalized() const // UnitQuaternions are always normalized, so checking for normalization is private to prevent abuse
            {
                static double const epsilon = 1e-14;
                double const ns = normSquared();
                return std::abs(1.0 - ns) < epsilon;
            }
        public:
            inline bool isValid() const
            {
                return isNormalized();
            }

            inline UnitQuaternion operator*(Vec3 const& otherV) const
            {
                REQUIRE(otherV.isNormalized());

                UnitQuaternion result(
                    /* m_r*0 */ -m_v.dot(otherV),
                    m_r*otherV + /* 0*m_v + */ m_v.cross(otherV)
                    );

                ENSURE(result.isValid());
                return result;
            }

            inline UnitQuaternion operator*(UnitQuaternion const& other) const
            {
                REQUIRE(other.isValid());

                UnitQuaternion result(
                    m_r*other.m_r - m_v.dot(other.m_v),
                    m_r*other.m_v + other.m_r*m_v + m_v.cross(other.m_v)
                    );

                ENSURE(result.isValid());
                return result;
            }

            inline UnitQuaternion conjugate() const
            {
                return UnitQuaternion(m_r, -m_v);
            }
            inline UnitQuaternion inverse() const // the inverse of a unit quaternion is its conjugate
            {
                return conjugate();
            }

            inline double getAngle() const
            {
                double const result = 2.0 * acos(m_r);
                return result;
            }

        private:
            inline Vec3 getAxis() const // this function is not currently used
            {
                CHECK(m_r > 0.0);
                REQUIRE(isValid());
                if (m_r == 1.0 || m_r == -1.0)
                    return Vec3(1.0, 0.0, 0.0); // any normalized axis will do
                else
                    return m_v / sqrt(1.0 - m_r*m_r);
            }
        public:

            //
            // because there is only a UnitQuaternion class, we can't tranform arbitatary vectors/points only
            // unit vectors/points on the unit sphere.
            //
            // TODO: Add regular Quaternion class and use it here in the UnitQuaternion class to support more than
            // just unit vector rotation
            //
            inline Vec3 rotateUnitVector(Vec3 const& p) const
            {
                REQUIRE(p.isNormalized());
                // rotation of vector p about quaternion q is p'=q*p*inverse(q). Note that p is promoted to a pure quaternion with r=0
                UnitQuaternion const q = *this * p * inverse();
                Vec3 result = q.m_v;
                ENSURE(result.isNormalized());
                return q.m_v;
            }
        };
    }

    namespace {
        /*
        * Calculate HeFF norm for one specific set of vectors (order and +/- direction.
        */
        double calculateSingleHeffNorm(
            double he1x, double he1y, double he1z,
            double he2x, double he2y, double he2z,
            double he3x, double he3y, double he3z,
            double ff1x, double ff1y, double ff1z,
            double ff2x, double ff2y, double ff2z,
            double ff3x, double ff3y, double ff3z)
        {
            double const epsilon = 1.0e-12;

            Vec3 const he1V = Vec3(he1x, he1y, he1z).normalize();
            Vec3 const he2V = Vec3(he2x, he2y, he2z).normalize();
            Vec3 const he3V = Vec3(he3x, he3y, he3z).normalize();

            Vec3 const ff1V = Vec3(ff1x, ff1y, ff1z).normalize();
            Vec3 const ff2V = Vec3(ff2x, ff2y, ff2z).normalize();
            Vec3 const ff3V = Vec3(ff3x, ff3y, ff3z).normalize();

            double norm = 0.0;
            // if the first vectors are not normalized now, we can assume they are zero vectors and just skip everything
            if (he1V.isNormalized() && ff1V.isNormalized())
            {
                // first rotate he1V to align with ff1V, but they might already be aligned
                Vec3 rotationAxis1;
                double angle1;
                double const he1DotFF1 = he1V.dot(ff1V);
                if (he1DotFF1 >= 1.0 - epsilon) // parallel-in-same-direction
                {
                    rotationAxis1 = Vec3(1.0, 0.0, 0.0); // any axis will do
                    angle1 = 0.0;
                }
                else if (he1DotFF1 <= -1.0 + epsilon) // parallel-in-opposite-direction
                {
                    // use any axis perpendicular to ff1V & he1V. Find one by crossing with either X or Y axis
                    if (std::abs(he1V.x()) < 0.5)
                        rotationAxis1 = Vec3(1.0, 0.0, 0.0).cross(ff1V).normalize();
                    else
                        rotationAxis1 = Vec3(0.0, 1.0, 0.0).cross(ff1V).normalize();
                    angle1 = M_PI;
                }
                else // needs rotation
                {
                    rotationAxis1 = he1V.cross(ff1V).normalize();
                    angle1 = acos(he1DotFF1);
                }
                UnitQuaternion r1Q(rotationAxis1, angle1); // TODO: look into creating the quaternion with the dot product directly to cos(angle/2)&sin(angle/2)

                // he1Q now rotates he1V into ff1V
                Vec3 const he1Vr = r1Q.rotateUnitVector(he1V); // he1V rotated
                CHECK(he1Vr.nearlyEquals(ff1V));
                // he1Q.inverse rotates ff1V into he1V
                UnitQuaternion const ir1Q = r1Q.inverse(); // inverse rotation 1 quaternion
                Vec3 const ff1Vir = ir1Q.rotateUnitVector(ff1V); // ff1V inverse rotated
                CHECK(ff1Vir.nearlyEquals(he1V));

                if (he2V.isNormalized() && ff2V.isNormalized())
                {
                    // rotate rotated-he2V to align with ff2V
                    Vec3 const he2rV = r1Q.rotateUnitVector(he2V); // get he2V-rotated

                    // rotated-he2V may already be aligned with ff2V
                    double angle2;
                    double he2rDotFF2 = he2rV.dot(ff2V);
                    if (he2rDotFF2 >= 1.0 - epsilon) // parallel-in-same-direction
                    {
                        angle2 = 0.0;
                    }
                    else if (he2rDotFF2 <= -1.0 + epsilon) // parallel-in-opposite-direction
                    {
                        angle2 = M_PI;
                    }
                    else // needs rotation
                    {
                        // instead of using ff1V, calculate rotation axis to get direction of rotation (and check our work so far)
                        Vec3 const cross2 = he2rV.cross(ff2V).normalize();
                        double const c2DotFF1 = cross2.dot(ff1V);
                        if (c2DotFF1 > 0.0)
                        {
                            CHECK(cross2.nearlyEquals(ff1V));
                            angle2 = acos(he2rDotFF2);
                        }
                        else
                        {
                            CHECK(cross2.nearlyEquals(-ff1V));
                            angle2 = -acos(he2rDotFF2); // rotate negatively around ff1V
                        }
                    }
                    UnitQuaternion r2Q(ff1V, angle2); // TODO: look into creating the quaternion with the dot product directly to cos(angle/2)&sin(angle/2)

                    // r2Q now rotates he2rV into ff2V
                    Vec3 const he2r2V = r2Q.rotateUnitVector(he2rV); // he2V rotated twice
                    CHECK(he2r2V.nearlyEquals(ff2V));
                    // r2Q.inverse rotates ff2V into he2rV
                    UnitQuaternion const ir2Q = r2Q.inverse(); // inverse rotation 2 quaternion
                    Vec3 const ff2irV = ir2Q.rotateUnitVector(ff2V); // ff2V inverse rotated
                    CHECK(ff2irV.nearlyEquals(he2rV));

                    // combine the two rotations
                    UnitQuaternion const combinedQ = r2Q * r1Q;

                    // check our work
                    Vec3 const cQhe1V = combinedQ.rotateUnitVector(he1V);
                    Vec3 const cQhe2V = combinedQ.rotateUnitVector(he2V);
                    double const d1 = cQhe1V.dot(ff1V);
                    double const d2 = cQhe2V.dot(ff2V);
                    CHECK(cQhe1V.nearlyEquals(ff1V));
                    CHECK(cQhe2V.nearlyEquals(ff2V));

                    // check third axis
                    if (he3V.isNormalized() && ff3V.isNormalized())
                    {
                        Vec3 const he3Vrc = combinedQ.rotateUnitVector(he3V); // he3V rotated by combined matrix

                        double const he3rcDotFF3 = he3Vrc.dot(ff3V);
                        if (he3rcDotFF3 > 0.0)
                        {
                            CHECK(he3Vrc.nearlyEquals(ff3V)); // parallel-in-same-direction
                            norm = combinedQ.getAngle();
                        }
                        else
                        {
                            CHECK(he3Vrc.nearlyEquals(-ff3V)); // parallel-in-opposite-direction
                            norm = M_PI; // maximum rotation angle possible
                        }
                    }
                    else
                    {
                        CHECK(he3V.nearlyZero() || ff3V.nearlyZero());
                        norm = combinedQ.getAngle();
                    }
                }
                else
                {
                    // Possible extension : when second vectors are zero, align third vectors
                    CHECK(he2V.nearlyZero() || ff2V.nearlyZero());
                    norm = angle1;
                }
            }
            else
            {
                CHECK(he1V.nearlyZero() || ff1V.nearlyZero());
                // Possible extension: when first vectors are zero, align second vectors
                norm = 0.0;
            }

            ENSURE(norm >= 0.0);
            return norm;
        }
    }

namespace {
    /*
    * Calculate HeFF norm for one orientation (vector order), but for any combination of +/- directions of the Hessian eigenvectors.
    */
    double calculateAnyDirectionHeffNorm(
        double he1x, double he1y, double he1z, // Hessian eigenvector 1
        double he2x, double he2y, double he2z, // Hessian eigenvector 2
        double he3x, double he3y, double he3z, // Hessian eigenvector 3
        double ff1x, double ff1y, double ff1z, // Frenet-Frame vector 1
        double ff2x, double ff2y, double ff2z, // Frenet-Frame vector 2
        double ff3x, double ff3y, double ff3z) // Frenet-Frame vector 3
    {
        double norm = 1000.0; // norm is an angle of rotation, so the maximum norm is pi

        for (int dir1 = -1; dir1 <= 1; dir1 += 2)
        {
            for (int dir2 = -1; dir2 <= 1; dir2 += 2)
            {
                for (int dir3 = -1; dir3 <= 1; dir3 += 2)
                {
                    double const dirNorm = calculateSingleHeffNorm(
                        dir1*he1x, dir1*he1y, dir1*he1z,
                        dir2*he2x, dir2*he2y, dir2*he2z,
                        dir3*he3x, dir3*he3y, dir3*he3z,
                        ff1x, ff1y, ff1z,
                        ff2x, ff2y, ff2z,
                        ff3x, ff3y, ff3z);
                    norm = std::min(dirNorm, norm);
                }
            }
        }
        return norm;
    }
}


/*
* Calculate the "norm" from the Hessian Eigensystem to the Frenet-Frame. This norm is the minimum
* rotation angle (in radians) required to rotate the Frenet into alignment with the Eigensystem, including
* all variations of that system using vectors in the opposite direction.
* The eigenvectors can be specified in any order, and the same for the Frenet-Frame vectors
* (The parenthetical comments below are just suggestions).
*/
double calculateHeffNorm(
    double he1x, double he1y, double he1z, // Hessian eigenvector 1
    double he2x, double he2y, double he2z, // Hessian eigenvector 2
    double he3x, double he3y, double he3z, // Hessian eigenvector 3
    double ff1x, double ff1y, double ff1z, // Frenet-Frame vector 1 (tangent)
    double ff2x, double ff2y, double ff2z, // Frenet-Frame vector 2 (normal)
    double ff3x, double ff3y, double ff3z) // Frenet-Frame vector 3 (binormal)
{
    double norm1 = calculateAnyDirectionHeffNorm(
        he1x, he1y, he1z,
        he2x, he2y, he2z,
        he3x, he3y, he3z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double norm2 = calculateAnyDirectionHeffNorm(
        he1x, he1y, he1z,
        he3x, he3y, he3z,
        he2x, he2y, he2z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double norm3 = calculateAnyDirectionHeffNorm(
        he2x, he2y, he2z,
        he1x, he1y, he1z,
        he3x, he3y, he3z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double norm4 = calculateAnyDirectionHeffNorm(
        he2x, he2y, he2z,
        he3x, he3y, he3z,
        he1x, he1y, he1z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double norm5 = calculateAnyDirectionHeffNorm(
        he3x, he3y, he3z,
        he1x, he1y, he1z,
        he2x, he2y, he2z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double norm6 = calculateAnyDirectionHeffNorm(
        he3x, he3y, he3z,
        he2x, he2y, he2z,
        he1x, he1y, he1z,
        ff1x, ff1y, ff1z,
        ff2x, ff2y, ff2z,
        ff3x, ff3y, ff3z);
    double const norm12 = std::min(norm1, norm2);
    double const norm34 = std::min(norm3, norm4);
    double const norm1234 = std::min(norm12, norm34);
    double const norm56 = std::min(norm5, norm6);
    double const norm1to6 = std::min(norm1234, norm56);

    return norm1to6;
}

}