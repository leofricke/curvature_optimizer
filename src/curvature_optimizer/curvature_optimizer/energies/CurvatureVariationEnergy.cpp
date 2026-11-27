#include "CurvatureVariationEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/DiffRules.h"
#include "../quantities/ScalarFieldGradientDerivatives.h"

CurvatureVariationEnergy::CurvatureVariationEnergy(CurvatureOptimizerWorkspace& ws, double pNorm)
    : GaussNewtonEnergyModule("CurvatureVariationEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_pNorm(pNorm)
{
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));
        if (!m_ws.m_mesh.is_boundary(v0) && !m_ws.m_mesh.is_boundary(v1) && !m_ws.m_mesh.is_boundary(v2))
        {
            m_faces.push_back(f);
        }
    }
}

static double absPowN(double x, double n)
{
    if (n == 1)
    {
        return std::abs(x);
    }
    const double a = std::max(1e-8, std::abs(x));
    return std::pow(a, n);
}

static double absPowNGradient(double x, double n, double& dx)
{
    const double sign = x >= 0.0 ? 1.0 : -1.0;
    if (n == 1)
    {
        dx = sign;
        return std::abs(x);
    }
    const double a = std::max(1e-8, std::abs(x));
    dx = (a <= 1e-8) ? 0.0 : sign * n * std::pow(a, n - 1.0);
    return std::pow(a, n);
}

size_t CurvatureVariationEnergy::numElements() const
{
    return m_faces.size();
}

size_t CurvatureVariationEnergy::numResiduals(size_t /*element*/) const
{
    return 6;
}

void CurvatureVariationEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const pmp::Face f = m_faces[element];
    const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
    const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));

    const IVertexCurvatureOperator& curv0 = m_ws.getCurrentVertexCurvature(v0);
    const IVertexCurvatureOperator& curv1 = m_ws.getCurrentVertexCurvature(v1);
    const IVertexCurvatureOperator& curv2 = m_ws.getCurrentVertexCurvature(v2);

    const double A = m_ws.m_currentFaceAreas[f];
    const double sqrtA = std::sqrt(A);

    if (!ctx.withJacobian())
    {
        const ScalarFieldFaceData fieldDataMax = { .m_face = f, .m_f0 = curv0.m_curv.m_kmax, .m_f1 = curv1.m_curv.m_kmax, .m_f2 = curv2.m_curv.m_kmax, .m_g0 = {}, .m_g1 = {}, .m_g2 = {} };
        const ScalarFieldFaceData fieldDataMin = { .m_face = f, .m_f0 = curv0.m_curv.m_kmin, .m_f1 = curv1.m_curv.m_kmin, .m_f2 = curv2.m_curv.m_kmin, .m_g0 = {}, .m_g1 = {}, .m_g2 = {} };
        const Eigen::Vector3d gradMax = computeScalarFieldGradient(m_ws, fieldDataMax);
        const Eigen::Vector3d gradMin = computeScalarFieldGradient(m_ws, fieldDataMin);
        for (size_t i = 0; i < 3; ++i)
        {
            ctx.residual(i) = sqrtA * absPowN(gradMax[i], 0.5 * m_pNorm);
            ctx.residual(3 + i) = sqrtA * absPowN(gradMin[i], 0.5 * m_pNorm);
        }
        return;
    }

    const ScalarFieldFaceData fieldDataMax = { .m_face = f,
                                               .m_f0 = curv0.m_curv.m_kmax,
                                               .m_f1 = curv1.m_curv.m_kmax,
                                               .m_f2 = curv2.m_curv.m_kmax,
                                               .m_g0 = m_ws.m_currentMaxCurvatureGradients[v0],
                                               .m_g1 = m_ws.m_currentMaxCurvatureGradients[v1],
                                               .m_g2 = m_ws.m_currentMaxCurvatureGradients[v2] };
    const ScalarFieldFaceData fieldDataMin = { .m_face = f,
                                               .m_f0 = curv0.m_curv.m_kmin,
                                               .m_f1 = curv1.m_curv.m_kmin,
                                               .m_f2 = curv2.m_curv.m_kmin,
                                               .m_g0 = m_ws.m_currentMinCurvatureGradients[v0],
                                               .m_g1 = m_ws.m_currentMinCurvatureGradients[v1],
                                               .m_g2 = m_ws.m_currentMinCurvatureGradients[v2] };

    LocalGradient<Eigen::Matrix3d, 20> Jkmax;
    LocalGradient<Eigen::Matrix3d, 20> Jkmin;
    const Eigen::Vector3d gradMax = computeScalarFieldGradientJacobian(m_ws, fieldDataMax, Jkmax);
    const Eigen::Vector3d gradMin = computeScalarFieldGradientJacobian(m_ws, fieldDataMin, Jkmin);

    LocalGradient<Eigen::RowVector3d> ga;
    applyChainRule1(A, [](double x) { return 0.5 / std::sqrt(x); }, m_ws.m_currentFaceAreaGradients[f], ga);

    // residuals 0-2: maximum curvature gradient, 3-5: minimum curvature gradient
    for (size_t part = 0; part < 2; ++part)
    {
        const Eigen::Vector3d& grad = part == 0 ? gradMax : gradMin;
        const LocalGradient<Eigen::Matrix3d, 20>& J = part == 0 ? Jkmax : Jkmin;
        for (size_t i = 0; i < 3; ++i)
        {
            double d;
            const double r = absPowNGradient(grad[i], 0.5 * m_pNorm, d);
            ctx.residual(3 * part + i) = sqrtA * r;

            JacobianRow& row = ctx.jacobianRow(3 * part + i);
            for (pmp::Vertex vv : J)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row.insertSegment(dofIdx * 3, sqrtA * d * J[vv].row(i));
                }
            }
            for (pmp::Vertex vv : ga)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row.insertSegment(dofIdx * 3, r * ga[vv]);
                }
            }
        }
    }
}
