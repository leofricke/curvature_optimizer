#include "CurvatureVariationEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/DiffRules.h"
#include "../quantities/ScalarFieldGradientDerivatives.h"

CurvatureVariationEnergy::CurvatureVariationEnergy(CurvatureOptimizerWorkspace& ws, double pNorm)
    : ModularGaussNewtonEnergyBlock(0, ws.m_dofs.numDofs()), m_ws(ws), m_pNorm(pNorm)
{
    m_name = "CurvatureVariationEnergy";

    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));
        if (m_ws.m_mesh.is_boundary(v0) || m_ws.m_mesh.is_boundary(v1) || m_ws.m_mesh.is_boundary(v2))
        {
            continue;
        }
        m_numResiduals += 6;
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

void CurvatureVariationEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));
        if (m_ws.m_mesh.is_boundary(v0) || m_ws.m_mesh.is_boundary(v1) || m_ws.m_mesh.is_boundary(v2))
        {
            continue;
        }

        IVertexCurvatureOperator& curv0 = m_ws.getCurrentVertexCurvature(v0);
        IVertexCurvatureOperator& curv1 = m_ws.getCurrentVertexCurvature(v1);
        IVertexCurvatureOperator& curv2 = m_ws.getCurrentVertexCurvature(v2);
        const ScalarFieldFaceData fieldDataMax = { .m_face = f, .m_f0 = curv0.m_curv.m_kmax, .m_f1 = curv1.m_curv.m_kmax, .m_f2 = curv2.m_curv.m_kmax, .m_g0 = {}, .m_g1 = {}, .m_g2 = {} };
        const ScalarFieldFaceData fieldDataMin = { .m_face = f, .m_f0 = curv0.m_curv.m_kmin, .m_f1 = curv1.m_curv.m_kmin, .m_f2 = curv2.m_curv.m_kmin, .m_g0 = {}, .m_g1 = {}, .m_g2 = {} };

        const Eigen::Vector3d gradMax = computeScalarFieldGradient(m_ws, fieldDataMax);
        const Eigen::Vector3d gradMin = computeScalarFieldGradient(m_ws, fieldDataMin);

        system.m_residual[ridx] = absPowN(gradMax[0], 0.5 * m_pNorm);
        system.m_residual[ridx + 1] = absPowN(gradMax[1], 0.5 * m_pNorm);
        system.m_residual[ridx + 2] = absPowN(gradMax[2], 0.5 * m_pNorm);
        system.m_residual[ridx + 3] = absPowN(gradMin[0], 0.5 * m_pNorm);
        system.m_residual[ridx + 4] = absPowN(gradMin[1], 0.5 * m_pNorm);
        system.m_residual[ridx + 5] = absPowN(gradMin[2], 0.5 * m_pNorm);

        const double A = m_ws.m_currentFaceAreas[f];
        const double sqrtA = std::sqrt(A);
        system.m_residual[ridx] *= sqrtA;
        system.m_residual[ridx + 1] *= sqrtA;
        system.m_residual[ridx + 2] *= sqrtA;
        system.m_residual[ridx + 3] *= sqrtA;
        system.m_residual[ridx + 4] *= sqrtA;
        system.m_residual[ridx + 5] *= sqrtA;

        ridx += 6;
    }
}

void CurvatureVariationEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row0 = system.createJacobianRow();
    spmat::SpVec<double> row1 = system.createJacobianRow();
    spmat::SpVec<double> row2 = system.createJacobianRow();
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));
        if (m_ws.m_mesh.is_boundary(v0) || m_ws.m_mesh.is_boundary(v1) || m_ws.m_mesh.is_boundary(v2))
        {
            continue;
        }

        IVertexCurvatureOperator& curv0 = m_ws.getCurrentVertexCurvature(v0);
        IVertexCurvatureOperator& curv1 = m_ws.getCurrentVertexCurvature(v1);
        IVertexCurvatureOperator& curv2 = m_ws.getCurrentVertexCurvature(v2);
        const LocalGradient<Eigen::RowVector3d, 8>& gkmax0 = m_ws.m_currentMaxCurvatureGradients[v0];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmin0 = m_ws.m_currentMinCurvatureGradients[v0];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmax1 = m_ws.m_currentMaxCurvatureGradients[v1];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmin1 = m_ws.m_currentMinCurvatureGradients[v1];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmax2 = m_ws.m_currentMaxCurvatureGradients[v2];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmin2 = m_ws.m_currentMinCurvatureGradients[v2];
        const ScalarFieldFaceData fieldDataMax = { .m_face = f, .m_f0 = curv0.m_curv.m_kmax, .m_f1 = curv1.m_curv.m_kmax, .m_f2 = curv2.m_curv.m_kmax, .m_g0 = gkmax0, .m_g1 = gkmax1, .m_g2 = gkmax2 };
        const ScalarFieldFaceData fieldDataMin = { .m_face = f, .m_f0 = curv0.m_curv.m_kmin, .m_f1 = curv1.m_curv.m_kmin, .m_f2 = curv2.m_curv.m_kmin, .m_g0 = gkmin0, .m_g1 = gkmin1, .m_g2 = gkmin2 };

        LocalGradient<Eigen::Matrix3d, 20> Jkmax;
        LocalGradient<Eigen::Matrix3d, 20> Jkmin;
        const Eigen::Vector3d gradMax = computeScalarFieldGradientJacobian(m_ws, fieldDataMax, Jkmax);
        const Eigen::Vector3d gradMin = computeScalarFieldGradientJacobian(m_ws, fieldDataMin, Jkmin);

        Eigen::RowVector3d dmax, dmin;
        system.m_residual[ridx] = absPowNGradient(gradMax[0], 0.5 * m_pNorm, dmax[0]);
        system.m_residual[ridx + 1] = absPowNGradient(gradMax[1], 0.5 * m_pNorm, dmax[1]);
        system.m_residual[ridx + 2] = absPowNGradient(gradMax[2], 0.5 * m_pNorm, dmax[2]);
        system.m_residual[ridx + 3] = absPowNGradient(gradMin[0], 0.5 * m_pNorm, dmin[0]);
        system.m_residual[ridx + 4] = absPowNGradient(gradMin[1], 0.5 * m_pNorm, dmin[1]);
        system.m_residual[ridx + 5] = absPowNGradient(gradMin[2], 0.5 * m_pNorm, dmin[2]);

        const double A = m_ws.m_currentFaceAreas[f];
        const double sqrtA = std::sqrt(A);
        LocalGradient<Eigen::RowVector3d> ga;
        applyChainRule1(A, [](double x) { return 0.5 / std::sqrt(x); }, m_ws.m_currentFaceAreaGradients[f], ga);

        const double r0 = system.m_residual[ridx];
        const double r1 = system.m_residual[ridx + 1];
        const double r2 = system.m_residual[ridx + 2];
        const double r3 = system.m_residual[ridx + 3];
        const double r4 = system.m_residual[ridx + 4];
        const double r5 = system.m_residual[ridx + 5];
        system.m_residual[ridx] *= sqrtA;
        system.m_residual[ridx + 1] *= sqrtA;
        system.m_residual[ridx + 2] *= sqrtA;
        system.m_residual[ridx + 3] *= sqrtA;
        system.m_residual[ridx + 4] *= sqrtA;
        system.m_residual[ridx + 5] *= sqrtA;

        {
            row0.clearPattern();
            row0.beginFill();
            row1.clearPattern();
            row1.beginFill();
            row2.clearPattern();
            row2.beginFill();
            for (pmp::Vertex vv : Jkmax)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row0.insertSegment(dofIdx * 3, sqrtA * dmax[0] * Jkmax[vv].row(0));
                    row1.insertSegment(dofIdx * 3, sqrtA * dmax[1] * Jkmax[vv].row(1));
                    row2.insertSegment(dofIdx * 3, sqrtA * dmax[2] * Jkmax[vv].row(2));
                }
            }
            for (pmp::Vertex vv : ga)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row0.insertSegment(dofIdx * 3, r0 * ga[vv]);
                    row1.insertSegment(dofIdx * 3, r1 * ga[vv]);
                    row2.insertSegment(dofIdx * 3, r2 * ga[vv]);
                }
            }
            row0.endFill();
            system.insertJacobianRow(ridx, row0);
            row1.endFill();
            system.insertJacobianRow(ridx + 1, row1);
            row2.endFill();
            system.insertJacobianRow(ridx + 2, row2);
        }

        {
            row0.clearPattern();
            row0.beginFill();
            row1.clearPattern();
            row1.beginFill();
            row2.clearPattern();
            row2.beginFill();
            for (pmp::Vertex vv : Jkmin)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row0.insertSegment(dofIdx * 3, sqrtA * dmin[0] * Jkmin[vv].row(0));
                    row1.insertSegment(dofIdx * 3, sqrtA * dmin[1] * Jkmin[vv].row(1));
                    row2.insertSegment(dofIdx * 3, sqrtA * dmin[2] * Jkmin[vv].row(2));
                }
            }
            for (pmp::Vertex vv : ga)
            {
                const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                if (dofIdx >= 0)
                {
                    row0.insertSegment(dofIdx * 3, r3 * ga[vv]);
                    row1.insertSegment(dofIdx * 3, r4 * ga[vv]);
                    row2.insertSegment(dofIdx * 3, r5 * ga[vv]);
                }
            }
            row0.endFill();
            system.insertJacobianRow(ridx + 3, row0);
            row1.endFill();
            system.insertJacobianRow(ridx + 4, row1);
            row2.endFill();
            system.insertJacobianRow(ridx + 5, row2);
        }

        ridx += 6;
    }
}

void CurvatureVariationEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const pmp::Halfedge h = m_ws.m_mesh.halfedge(f);
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(h);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(h);
        const pmp::Vertex v2 = m_ws.m_mesh.to_vertex(m_ws.m_mesh.next_halfedge(h));
        if (m_ws.m_mesh.is_boundary(v0) || m_ws.m_mesh.is_boundary(v1) || m_ws.m_mesh.is_boundary(v2))
        {
            continue;
        }

        weights.segment<6>(ridx).setConstant(energyWeight);
        ridx += 6;
    }
}