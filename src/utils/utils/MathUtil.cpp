#include "MathUtil.h"

pmp::dmat2 rotationMatrix(double angle)
{
    const double c = std::cos(angle);
    const double s = std::sin(angle);

    pmp::dmat2 m(0.0);
    m(0, 0) = c;
    m(0, 1) = -s;
    m(1, 0) = s;
    m(1, 1) = c;
    return m;
}

pmp::dmat4 rotationMatrix(const pmp::dvec3& axis, double angle)
{
    pmp::dmat4 m(0.0);
    const double a = angle;
    const double c = std::cos(a);
    const double s = std::sin(a);
    const double one_m_c = 1.0 - c;
    const pmp::dvec3 ax = normalize(axis);

    m(0, 0) = ax[0] * ax[0] * one_m_c + c;
    m(0, 1) = ax[0] * ax[1] * one_m_c - ax[2] * s;
    m(0, 2) = ax[0] * ax[2] * one_m_c + ax[1] * s;

    m(1, 0) = ax[1] * ax[0] * one_m_c + ax[2] * s;
    m(1, 1) = ax[1] * ax[1] * one_m_c + c;
    m(1, 2) = ax[1] * ax[2] * one_m_c - ax[0] * s;

    m(2, 0) = ax[2] * ax[0] * one_m_c - ax[1] * s;
    m(2, 1) = ax[2] * ax[1] * one_m_c + ax[0] * s;
    m(2, 2) = ax[2] * ax[2] * one_m_c + c;

    m(3, 3) = 1.0;

    return m;
}

pmp::dvec3 affineTransform(const pmp::dmat4& m, const pmp::dvec3& v)
{
    const double x = m(0, 0) * v[0] + m(0, 1) * v[1] + m(0, 2) * v[2] + m(0, 3);
    const double y = m(1, 0) * v[0] + m(1, 1) * v[1] + m(1, 2) * v[2] + m(1, 3);
    const double z = m(2, 0) * v[0] + m(2, 1) * v[1] + m(2, 2) * v[2] + m(2, 3);
    return pmp::dvec3(x, y, z);
}

pmp::dvec3 linearTransform(const pmp::dmat4& m, const pmp::dvec3& v)
{
    const double x = m(0, 0) * v[0] + m(0, 1) * v[1] + m(0, 2) * v[2];
    const double y = m(1, 0) * v[0] + m(1, 1) * v[1] + m(1, 2) * v[2];
    const double z = m(2, 0) * v[0] + m(2, 1) * v[1] + m(2, 2) * v[2];
    return pmp::dvec3(x, y, z);
}

pmp::dmat2 inverse(const pmp::dmat2& m)
{
    const double det = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
    if (std::abs(det) < 1e-20)
    {
        return pmp::dmat2(std::numeric_limits<double>::infinity());
    }

    pmp::dmat2 inv;
    inv(0, 0) = m(1, 1);
    inv(0, 1) = -m(0, 1);
    inv(1, 0) = -m(1, 0);
    inv(1, 1) = m(0, 0);
    inv /= det;
    return inv;
}
