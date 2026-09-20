#include "GaussianCurvatureDerivatives.h"

#include "DiffRules.h"
#include "../CurvatureOptimizerWorkspace.h"
#include <numbers>

double computeVertexIntegratedGaussianCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    if (ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    double K = 2.0 * std::numbers::pi;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        K -= ws.m_currentInteriorAngles[h];
    }

    return K;
}

double computeVertexGaussianCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    if (ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    return computeVertexIntegratedGaussianCurvature(ws, v) / ws.m_currentVertexDualAreas[v];
}

double computeVertexIntegratedGaussianCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g)
{
    if (ws.m_mesh.is_boundary(v))
    {
        return 0.0;
    }

    double K = 2.0 * std::numbers::pi;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            K -= ws.m_currentInteriorAngles[h];

            const LocalGradient<Eigen::RowVector3d, 3>& gb = ws.m_currentInteriorAngleGradients[h];

            g -= gb;
        }
    }

    return K;
}

double computeVertexGaussianCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g)
{
    const double K = computeVertexIntegratedGaussianCurvatureGradient(ws, v, g);
    if (ws.m_mesh.is_boundary(v))
    {
        return K;
    }

    const double A = ws.m_currentVertexDualAreas[v];
    const LocalGradient<Eigen::RowVector3d, 8>& ga = ws.m_currentVertexDualAreaGradients[v];

    applyQuotientRule1Inplace(K, A, ga, g);

    return K / A;
}
