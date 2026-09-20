#include "BoundaryCurvatureDerivatives.h"

#include "DiffRules.h"
#include "../CurvatureOptimizerWorkspace.h"
#include <numbers>

double computeIntegratedBoundaryCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    if (!ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    double k = std::numbers::pi;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        k -= ws.m_currentInteriorAngles[h];
    }

    return k;
}

double computeBoundaryCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    if (!ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    double L = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (ws.m_mesh.is_boundary(ws.m_mesh.edge(h)))
        {
            const Eigen::Vector3d p0 = ws.position(ws.m_mesh.from_vertex(h));
            const Eigen::Vector3d p1 = ws.position(ws.m_mesh.to_vertex(h));
            const Eigen::Vector3d e = p1 - p0;
            const double le = e.norm();

            L += 0.5 * le;
        }
    }

    return computeIntegratedBoundaryCurvature(ws, v) / L;
}

double computeIntegratedBoundaryCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d>& g)
{
    if (!ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    double k = std::numbers::pi;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            k -= ws.m_currentInteriorAngles[h];

            const LocalGradient<Eigen::RowVector3d, 3>& gb = ws.m_currentInteriorAngleGradients[h];

            for (pmp::Vertex vv : gb)
            {
                g[vv] -= gb[vv];
            }
        }
    }

    return k;
}

double computeBoundaryCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d>& g)
{
    const double k = computeIntegratedBoundaryCurvatureGradient(ws, v, g);
    if (!ws.m_mesh.is_boundary(v))
    {
        return k;
    }

    double L = 0.0;
    LocalGradient<Eigen::RowVector3d> gl;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (ws.m_mesh.is_boundary(ws.m_mesh.edge(h)))
        {
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v);
            const Eigen::Vector3d p1 = ws.position(v1);
            const Eigen::Vector3d e = p1 - p0;
            const double le = e.norm();
            const Eigen::Vector3d e0 = e.normalized();

            L += 0.5 * le;

            gl[v] -= 0.5 * e0;
            gl[v1] += 0.5 * e0;
        }
    }

    applyQuotientRule1Inplace(k, L, gl, g);

    return k / L;
}