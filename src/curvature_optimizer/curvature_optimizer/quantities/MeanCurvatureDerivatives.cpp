#include "MeanCurvatureDerivatives.h"

#include "DiffRules.h"
#include "../CurvatureOptimizerWorkspace.h"

double computeVertexIntegratedMeanCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    double H = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(ws.m_mesh.edge(h)))
        {
            const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);

            const Eigen::Vector3d ev = p1 - p0;
            const double l = ev.norm();

            const double beta = ws.m_currentDihedralAngles[ws.m_mesh.edge(h)];
            H += beta * 0.25 * l;
        }
    }

    return H;
}

double computeVertexMeanCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    return computeVertexIntegratedMeanCurvature(ws, v) / ws.m_currentVertexDualAreas[v];
}

double computeVertexIntegratedMeanCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g)
{
    double H = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(ws.m_mesh.edge(h)))
        {
            const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);

            const Eigen::Vector3d ev = p1 - p0;
            const double l = ev.norm();
            const Eigen::Vector3d ev0 = ev.normalized();

            const double beta = ws.m_currentDihedralAngles[ws.m_mesh.edge(h)];
            H += beta * 0.25 * l;

            // dihedral angle term
            const LocalGradient<Eigen::RowVector3d, 4>& gb = ws.m_currentDihedralAngleGradients[ws.m_mesh.edge(h)];

            // length term
            LocalGradient<Eigen::RowVector3d> gl;
            gl[v1] = ev0;
            gl[v0] = -ev0;

            applyProductRule1(beta, l, gb, gl, g);
        }
    }
    g *= 0.25;

    return H;
}

double computeVertexMeanCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g)
{
    const double H = computeVertexIntegratedMeanCurvatureGradient(ws, v, g);

    const double A = ws.m_currentVertexDualAreas[v];
    const LocalGradient<Eigen::RowVector3d, 8>& ga = ws.m_currentVertexDualAreaGradients[v];

    applyQuotientRule1Inplace(H, A, ga, g);

    return H / A;
}
