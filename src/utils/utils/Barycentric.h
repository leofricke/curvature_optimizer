#pragma once

#include <pmp/surface_mesh.h>

struct BarycentricCoords
{
    BarycentricCoords()
        : m_alpha(0.0), m_beta(0.0), m_gamma(0.0)
    {
    }

    BarycentricCoords(double alpha, double beta)
        : m_alpha(alpha), m_beta(beta), m_gamma(1.0 - m_alpha - m_beta)
    {
    }

    BarycentricCoords(double alpha, double beta, double gamma)
        : m_alpha(alpha), m_beta(beta), m_gamma(gamma)
    {
    }

    bool isValid() const
    {
        return std::abs(m_alpha + m_beta + m_gamma - 1.0) < 1e-8;
    }

    bool isInside() const
    {
        return m_alpha >= 0.0 && m_beta >= 0.0 && m_gamma >= 0.0;
    }

    template <class T>
    T interpolate(const T& a, const T& b, const T& c) const
    {
        return m_alpha * a + m_beta * b + m_gamma * c;
    }

    BarycentricCoords& operator+=(const BarycentricCoords& bary)
    {
        m_alpha += bary.m_alpha;
        m_beta += bary.m_beta;
        m_gamma += bary.m_gamma;
        return *this;
    }

    BarycentricCoords& operator-=(const BarycentricCoords& bary)
    {
        m_alpha -= bary.m_alpha;
        m_beta -= bary.m_beta;
        m_gamma -= bary.m_gamma;
        return *this;
    }

    double m_alpha;
    double m_beta;
    double m_gamma;
};

inline BarycentricCoords operator+(const BarycentricCoords& bary0, const BarycentricCoords& bary1)
{
    BarycentricCoords p = bary0;
    p += bary1;
    return p;
}

inline BarycentricCoords operator-(const BarycentricCoords& bary0, const BarycentricCoords& bary1)
{
    BarycentricCoords p = bary0;
    p -= bary1;
    return p;
}

struct BarycentricPoint : public BarycentricCoords
{
    BarycentricPoint() = default;

    BarycentricPoint(pmp::Halfedge h, double alpha, double beta)
        : BarycentricCoords(alpha, beta), m_halfedge(h)
    {
    }

    BarycentricPoint(pmp::Halfedge h, double alpha, double beta, double gamma)
        : BarycentricCoords(alpha, beta, gamma), m_halfedge(h)
    {
    }

    explicit BarycentricPoint(pmp::Vertex v, const pmp::SurfaceMesh& mesh)
        : BarycentricCoords(1.0, 0.0, 0.0)
    {
        for (pmp::Halfedge h : mesh.halfedges(v))
        {
            if (!mesh.is_boundary(h))
            {
                m_halfedge = h;
                break;
            }
        }
    }

    pmp::Vertex aVertex(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.from_vertex(m_halfedge);
    }

    pmp::Vertex bVertex(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.to_vertex(m_halfedge);
    }

    pmp::Vertex cVertex(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.to_vertex(mesh.next_halfedge(m_halfedge));
    }

    pmp::Face face(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.face(m_halfedge);
    }

    const pmp::Point& aPoint(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.position(aVertex(mesh));
    }

    const pmp::Point& bPoint(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.position(bVertex(mesh));
    }

    const pmp::Point& cPoint(const pmp::SurfaceMesh& mesh) const
    {
        return mesh.position(cVertex(mesh));
    }

    pmp::Point point(const pmp::SurfaceMesh& mesh) const
    {
        const pmp::Point& a = mesh.position(aVertex(mesh));
        const pmp::Point& b = mesh.position(bVertex(mesh));
        const pmp::Point& c = mesh.position(cVertex(mesh));

        return interpolate(a, b, c);
    }

    pmp::Normal normal(const pmp::SurfaceMesh& mesh) const
    {
        const auto normals = mesh.get_vertex_property<pmp::Normal>("v:normal");
        const pmp::Normal& na = normals[aVertex(mesh)];
        const pmp::Normal& nb = normals[bVertex(mesh)];
        const pmp::Normal& nc = normals[cVertex(mesh)];

        return interpolate(na, nb, nc);
    }

    static BarycentricPoint compute(const pmp::SurfaceMesh& mesh, pmp::Halfedge h, const pmp::Point& p)
    {
#ifndef NDEBUG
        if (mesh.valence(mesh.face(h)) != 3)
        {
            throw std::runtime_error("BarycentricPoint only works on triangles");
        }
#endif

        const pmp::dvec3 a(mesh.position(mesh.from_vertex(h)));
        const pmp::dvec3 b(mesh.position(mesh.to_vertex(h)));
        const pmp::dvec3 c(mesh.position(mesh.to_vertex(mesh.next_halfedge(h))));

        const pmp::dvec3 va(a - c);
        const pmp::dvec3 vb(b - c);
        const pmp::dvec3 vp(static_cast<pmp::dvec3>(p) - c);

        const double d00 = dot(va, va);
        const double d01 = dot(va, vb);
        const double d11 = dot(vb, vb);
        const double d20 = dot(vp, va);
        const double d21 = dot(vp, vb);

        const double denom = d00 * d11 - d01 * d01;

        const double alpha = (d11 * d20 - d01 * d21) / denom;
        const double beta = (d00 * d21 - d01 * d20) / denom;
        const double gamma = 1.0 - alpha - beta;

        return BarycentricPoint(h, alpha, beta, gamma);
    }

    BarycentricPoint& operator+=(const BarycentricPoint& bary)
    {
        m_alpha += bary.m_alpha;
        m_beta += bary.m_beta;
        m_gamma += bary.m_gamma;
        return *this;
    }

    BarycentricPoint& operator-=(const BarycentricPoint& bary)
    {
        m_alpha -= bary.m_alpha;
        m_beta -= bary.m_beta;
        m_gamma -= bary.m_gamma;
        return *this;
    }

    pmp::Halfedge m_halfedge;
};

inline BarycentricPoint operator+(const BarycentricPoint& bary0, const BarycentricPoint& bary1)
{
    BarycentricPoint p = bary0;
    p += bary1;
    return p;
}

inline BarycentricPoint operator-(const BarycentricPoint& bary0, const BarycentricPoint& bary1)
{
    BarycentricPoint p = bary0;
    p -= bary1;
    return p;
}