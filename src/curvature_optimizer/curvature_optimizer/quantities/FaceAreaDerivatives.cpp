#include "FaceAreaDerivatives.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "utils/EigenMathUtil.h"

double computeFaceArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    if (ws.m_mesh.valence(f) == 3)
    {
        return computeTriangleArea(ws, f);
    }
    else
    {
        return computePolygonArea(ws, f);
    }
}

double computeFaceAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g)
{
    if (ws.m_mesh.valence(f) == 3)
    {
        return computeTriangleAreaGradient(ws, f, g);
    }
    else
    {
        return computePolygonAreaGradient(ws, f, g);
    }
}

double computeTriangleArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    const pmp::Halfedge h = ws.m_mesh.halfedge(f);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;

    return 0.5 * e01.cross(-e20).norm();
}

double computeTriangleAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g)
{
    const pmp::Halfedge h = ws.m_mesh.halfedge(f);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e12 = p2 - p1;
    const Eigen::Vector3d e20 = p0 - p2;

    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();
    g[v0] = 0.5 * n.cross(e12);
    g[v1] = 0.5 * n.cross(e20);
    g[v2] = 0.5 * n.cross(e01);

    return 0.5 * doubleArea;
}

double computePolygonArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p1 = ws.position(v1);
        n += p0.cross(p1);
    }

    return 0.5 * n.norm();
}

double computePolygonAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g)
{
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p1 = ws.position(v1);
        n += p0.cross(p1);
    }
    const double doubleArea = n.norm();
    n.normalize();

    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p2 = ws.position(v2);

        g[v1] = n.cross(p0 - p2);
    }
    g *= 0.5;

    return 0.5 * doubleArea;
}

Eigen::Vector3d computeFaceNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    if (ws.m_mesh.valence(f) == 3)
    {
        return computeTriangleNormal(ws, f);
    }
    else
    {
        return computePolygonNormal(ws, f);
    }
}

Eigen::Vector3d computeFaceNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d>& J)
{
    if (ws.m_mesh.valence(f) == 3)
    {
        return computeTriangleNormalJacobian(ws, f, J);
    }
    else
    {
        return computePolygonNormalJacobian(ws, f, J);
    }
}

Eigen::Vector3d computeTriangleNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    const pmp::Halfedge h = ws.m_mesh.halfedge(f);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;

    Eigen::Vector3d n = e01.cross(-e20);
    n.normalize();

    return n;
}

Eigen::Vector3d computeTriangleNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d>& J)
{
    const pmp::Halfedge h = ws.m_mesh.halfedge(f);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e12 = p2 - p1;
    const Eigen::Vector3d e20 = p0 - p2;

    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();

    const Eigen::Matrix3d P = (Eigen::Matrix3d::Identity() - n * n.transpose()) / doubleArea;
    J[v0] = P * crossMatrix<double>(e12);
    J[v1] = P * crossMatrix<double>(e20);
    J[v2] = P * crossMatrix<double>(e01);

    return n;
}

Eigen::Vector3d computePolygonNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p1 = ws.position(v1);
        n += p0.cross(p1);
    }
    n.normalize();

    return n;
}

Eigen::Vector3d computePolygonNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d>& J)
{
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p1 = ws.position(v1);
        n += p0.cross(p1);
    }
    const double doubleArea = n.norm();
    n.normalize();

    for (pmp::Halfedge h : ws.m_mesh.halfedges(f))
    {
        const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
        const Eigen::Vector3d p0 = ws.position(v0);
        const Eigen::Vector3d p2 = ws.position(v2);

        J[v1] = crossMatrix<double>(p0 - p2);
    }

    const Eigen::Matrix3d P = (Eigen::Matrix3d::Identity() - n * n.transpose()) / doubleArea;
    for (pmp::Vertex v : J)
    {
        J[v] = P * J[v];
    }

    return n;
}
