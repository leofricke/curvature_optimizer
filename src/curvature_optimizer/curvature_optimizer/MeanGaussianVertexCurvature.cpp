#include "MeanGaussianVertexCurvature.h"

#include "CurvatureOptimizerWorkspace.h"
#include "quantities/DiffRules.h"
#include "quantities/GaussianCurvatureDerivatives.h"
#include "quantities/MeanCurvatureDerivatives.h"
#include "utils/RobustSqrt.h"

static SimpleRobustSqrt s_robustSqrt(RobustSqrt::Mode::Clamp, 1e-8);

MeanGaussianVertexCurvature::MeanGaussianVertexCurvature()
{
}

Eigen::Vector2d MeanGaussianVertexCurvature::computePrincipalCurvatures(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    double H = computeVertexIntegratedMeanCurvature(ws, v);
    double K = computeVertexIntegratedGaussianCurvature(ws, v);
    const double A = ws.m_currentVertexDualAreas[v];
    H /= A;
    K /= A;

    const double d = H * H - K;
    const double c = s_robustSqrt.evaluate(d);

    const double kmax = H + c;
    const double kmin = H - c;

    m_curv.m_kmax = kmax;
    m_curv.m_kmin = kmin;
    return { m_curv.m_kmax, m_curv.m_kmin };
}

Eigen::Vector2d MeanGaussianVertexCurvature::computePrincipalCurvaturesGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& gkmax,
                                                                                LocalGradient<Eigen::RowVector3d, 8>& gkmin)
{
    LocalGradient<Eigen::RowVector3d, 8> gh;
    double H = computeVertexIntegratedMeanCurvatureGradient(ws, v, gh);

    LocalGradient<Eigen::RowVector3d, 8> gk;
    double K = computeVertexIntegratedGaussianCurvatureGradient(ws, v, gk);

    const double A = ws.m_currentVertexDualAreas[v];
    const LocalGradient<Eigen::RowVector3d, 8>& ga = ws.m_currentVertexDualAreaGradients[v];
    applyQuotientRule1Inplace(H, A, ga, gh);
    applyQuotientRule1Inplace(K, A, ga, gk);
    H /= A;
    K /= A;

    const double d = H * H - K;
    const double c = s_robustSqrt.evaluate(d);

    gkmax += gh;
    gkmin += gh;

    LocalGradient<Eigen::RowVector3d, 8> gc;
    applyChainRule1(H, [](double x) { return 2.0 * x; }, gh, gc);

    gc -= gk;

    applyChainRule1Inplace(d, [&](double x) { return s_robustSqrt.evaluateDeriv1(x); }, gc);

    gkmax += gc;
    gkmin -= gc;

    const double kmax = H + c;
    const double kmin = H - c;

    m_curv.m_kmax = kmax;
    m_curv.m_kmin = kmin;
    return { m_curv.m_kmax, m_curv.m_kmin };
}
