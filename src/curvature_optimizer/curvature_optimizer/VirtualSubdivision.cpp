#include "VirtualSubdivision.h"

static Eigen::Matrix<double, Eigen::Dynamic, 3> polygonMatrix(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    const size_t val = ws.m_mesh.valence(f);
    Eigen::Matrix<double, Eigen::Dynamic, 3> poly(val, 3);
    size_t i = 0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        poly.row(i) = ws.position(ws.m_mesh.from_vertex(h));
        i++;
    }
    return poly;
}

void VirtualSubdivision::updateCenterWeights(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    m_centerWeights.clear();
    const Eigen::Matrix<double, Eigen::Dynamic, 3> poly = polygonMatrix(ws, f);
    const int n = poly.rows();
    if (n == 3)
    {
        m_centerWeights.push_back(0.0);
        m_centerWeights.push_back(0.0);
        m_centerWeights.push_back(1.0);
        return;
    }

    // setup array of positions and edges
    std::vector<Eigen::Vector3d> d(n);
    for (int i = 0; i < n; ++i)
    {
        d[i] = poly.row((i + 1) % n) - poly.row(i);
    }

    // setup matrix A and rhs b
    // see Equation (38) of "Polygon Laplacian made simple", Eurographics 2020
    Eigen::MatrixXd A(n + 1, n);
    Eigen::VectorXd b(n + 1);
    for (int j = 0; j < n; ++j)
    {
        for (int i = j; i < n; ++i)
        {
            double Aij(0.0), bi(0.0);
            for (int k = 0; k < n; ++k)
            {
                Aij += poly.row(j).cross(d[k]).dot(poly.row(i).cross(d[k]));
                bi += poly.row(i).cross(d[k]).dot(poly.row(k).cross(d[k]));
            }
            A(i, j) = A(j, i) = Aij;
            b(i) = bi;
        }
    }
    for (int j = 0; j < n; ++j)
    {
        A(n, j) = 1.0;
    }
    b(n) = 1.0;

    const Eigen::VectorXd solution = A.completeOrthogonalDecomposition().solve(b).topRows(n);
    for (int i = 0; i < n; ++i)
    {
        m_centerWeights.push_back(solution[i]);
    }
}

Eigen::Vector3d VirtualSubdivision::getCenterPoint(const CurvatureOptimizerWorkspace& ws, pmp::Face f) const
{
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    size_t idx = 0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        center += m_centerWeights[idx] * ws.position(ws.m_mesh.from_vertex(h));
        idx++;
    }
    return center;
}

double VirtualSubdivision::getVertexWeight(const CurvatureOptimizerWorkspace& ws, pmp::Face f, pmp::Vertex v) const
{
    size_t hidx = 0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        if (ws.m_mesh.from_vertex(h) == v)
        {
            return m_centerWeights[hidx];
        }
        hidx++;
    }
    throw std::runtime_error("Vertex not found");
}

void VirtualSubdivision::initVirtualTriangles(const CurvatureOptimizerWorkspace& ws, pmp::Face f, InlineVector<VirtualTriangle, 4>& triangles) const
{
    triangles.clear();
    const size_t val = ws.m_mesh.valence(f);
    const Eigen::Vector3d center = getCenterPoint(ws, f);

    if (val == 3)
    {
        const pmp::Halfedge h = ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p1 = ws.position(v1);

        triangles.emplace_back(h, p0, p1, center);
    }
    else
    {
        for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
        {
            const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);

            triangles.emplace_back(h, p0, p1, center);
        }
    }
}

VirtualTriangle::VirtualTriangle(pmp::Halfedge h, const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& center)
{
    m_halfedge = h;

    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e12 = center - p1;
    const Eigen::Vector3d e20 = p0 - center;

    m_le01 = e01.norm();
    m_le12 = e12.norm();
    m_le20 = e20.norm();
    m_area = 0.5 * e01.cross(-e20).norm();

    const double a0 = e12.norm();
    const double a1 = e20.norm();
    const double a2 = e01.norm();
    const double s = a0 + a1 + a2;

    m_incenter = BarycentricCoords(a0 / s, a1 / s, a2 / s);

    const double ha1 = 0.5 * s - a2;
    const double ha2 = 0.5 * s - a1;
    const double has = ha1 + ha2;
    m_q12 = BarycentricCoords(0.0, ha1 / has, ha2 / has);

    const double hb0 = 0.5 * s - a2;
    const double hb2 = 0.5 * s - a0;
    const double hbs = hb0 + hb2;
    m_q20 = BarycentricCoords(hb0 / hbs, 0.0, hb2 / hbs);

    const double hc0 = 0.5 * s - a1;
    const double hc1 = 0.5 * s - a0;
    const double hcs = hc0 + hc1;
    m_q01 = BarycentricCoords(hc0 / hcs, hc1 / hcs, 0.0);
}
