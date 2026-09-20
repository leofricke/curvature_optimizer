#pragma once

#include <pmp/surface_mesh.h>

struct Ray
{
    Ray()
        : m_origin(0.0), m_dir(0.0)
    {
    }

    Ray(const pmp::vec3& origin, const pmp::vec3& dir)
        : m_origin(origin), m_dir(normalize(dir))
    {
    }

    pmp::vec3 evaluate(double t) const
    {
        return m_origin + t * m_dir;
    }

    pmp::Scalar distanceToPoint(const pmp::vec3& point, pmp::Scalar& t, pmp::vec3& rayPos) const
    {
        t = std::max(0.0f, dot(m_dir, point - m_origin));
        rayPos = evaluate(t);
        return distance(rayPos, point);
    }

    bool intersectTriangle(const pmp::vec3& a, const pmp::vec3& b, const pmp::vec3& c, pmp::Scalar& tHit, pmp::vec3& pointHit) const
    {
        const pmp::vec3 e0 = b - a;
        const pmp::vec3 e1 = c - a;

        const pmp::vec3 pvec = cross(m_dir, e1);
        const pmp::Scalar det = dot(e0, pvec);
        if (std::abs(det) < 1e-8f)
        {
            return false;
        }

        const pmp::Scalar invDet = 1.0f / det;
        const pmp::vec3 tvec = m_origin - a;
        const pmp::vec3 qvec = cross(tvec, e0);
        const pmp::Scalar u = dot(tvec, pvec) * invDet;
        const pmp::Scalar v = dot(m_dir, qvec) * invDet;
        tHit = dot(e1, qvec) * invDet;

        if (tHit >= 0.0f && u >= 0.0f && v >= 0.0f && u + v <= 1.0f)
        {
            pointHit = evaluate(tHit);
            return true;
        }

        return false;
    }

    bool intersectSphere(const pmp::vec3& center, pmp::Scalar radius, pmp::Scalar& tHit, pmp::vec3& pointHit) const
    {
        const pmp::vec3 oc = m_origin - center;
        const pmp::Scalar b = dot(m_dir, oc);
        const pmp::Scalar c = dot(oc, oc) - radius * radius;
        const pmp::Scalar discriminant = b * b - c;
        if (discriminant < 0.0)
        {
            return false;
        }
        const pmp::Scalar sqrtD = std::sqrt(discriminant);

        // two solutions
        pmp::Scalar t0 = -b - sqrtD;
        pmp::Scalar t1 = -b + sqrtD;
        if (t1 < 0.0)
        {
            return false;
        }

        tHit = (t0 >= 0.0) ? t0 : t1;
        pointHit = evaluate(tHit);

        return true;
    }

    pmp::vec3 m_origin;
    pmp::vec3 m_dir;
};
