#include "CurvatureEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"

CurvatureEnergy::CurvatureEnergy(CurvatureOptimizerWorkspace& ws, const std::vector<TargetVertexCurvature>& targetVertexCurvatures)
    : ModularGaussNewtonEnergyBlock(targetVertexCurvatures.size() * 2, ws.m_dofs.numDofs()), m_ws(ws), m_targetVertexCurvatures(targetVertexCurvatures)
{
    m_name = "CurvatureEnergy";
}

void CurvatureEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (const TargetVertexCurvature& target : m_targetVertexCurvatures)
    {
        IVertexCurvatureOperator& currentCurv = m_ws.getCurrentVertexCurvature(target.m_vertex);

        const Eigen::Vector2d targetCurv = { target.m_kmax, target.m_kmin };
        system.m_residual.segment<2>(ridx) = currentCurv.computePrincipalCurvatures(m_ws, target.m_vertex) - targetCurv;
        ridx += 2;
    }
}

void CurvatureEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row = system.createJacobianRow();
    for (const TargetVertexCurvature& target : m_targetVertexCurvatures)
    {
        const Eigen::Vector2d targetCurv = { target.m_kmax, target.m_kmin };
        IVertexCurvatureOperator& currentCurv = m_ws.getCurrentVertexCurvature(target.m_vertex);
        system.m_residual.segment<2>(ridx) = currentCurv.computePrincipalCurvatures(m_ws, target.m_vertex) - targetCurv;

        const LocalGradient<Eigen::RowVector3d, 8>& gkmax = m_ws.m_currentMaxCurvatureGradients[target.m_vertex];
        const LocalGradient<Eigen::RowVector3d, 8>& gkmin = m_ws.m_currentMinCurvatureGradients[target.m_vertex];

        row.clearPattern();
        row.beginFill();
        for (pmp::Vertex vv : gkmax)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.insertSegment(dofIdx * 3, gkmax[vv]);
            }
        }
        row.endFill();
        system.insertJacobianRow(ridx, row);

        row.clearPattern();
        row.beginFill();
        for (pmp::Vertex vv : gkmin)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.insertSegment(dofIdx * 3, gkmin[vv]);
            }
        }
        row.endFill();
        system.insertJacobianRow(ridx + 1, row);

        ridx += 2;
    }
}

void CurvatureEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (const TargetVertexCurvature& target : m_targetVertexCurvatures)
    {
        weights.segment<2>(ridx).setConstant(energyWeight * m_ws.m_currentVertexDualAreas[target.m_vertex]);
        ridx += 2;
    }
}
