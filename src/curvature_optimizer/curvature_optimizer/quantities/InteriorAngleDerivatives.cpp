#include "InteriorAngleDerivatives.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "utils/EigenMathUtil.h"

double computeInteriorAngle(const CurvatureOptimizerWorkspace& ws, pmp::Halfedge h)
{
    if (ws.m_mesh.is_boundary(h))
    {
        return 0.0;
    }
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.from_vertex(ws.m_mesh.prev_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;
    const double le01 = e01.squaredNorm();
    const double le20 = e20.squaredNorm();
    const double d = e01.dot(-e20) / std::sqrt(le01 * le20);

    return std::acos(std::clamp(d, -1.0, 1.0));
}

double computeInteriorAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Halfedge h, LocalGradient<Eigen::RowVector3d, 3>& g)
{
    if (ws.m_mesh.is_boundary(h))
    {
        return 0.0;
    }
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.from_vertex(ws.m_mesh.prev_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;
    Eigen::Vector3d n = e01.cross(-e20);
    n.normalize();
    const double le01 = e01.squaredNorm();
    const double le20 = e20.squaredNorm();
    const double d = e01.dot(-e20) / std::sqrt(le01 * le20);

    g[v1] = e01.cross(n) / le01;
    g[v2] = e20.cross(n) / le20;
    g[v0] = -g[v1] - g[v2];

    return std::acos(std::clamp(d, -1.0, 1.0));
}
