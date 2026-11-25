#include "BoundaryCurvatureEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "../quantities/BoundaryCurvatureDerivatives.h"

BoundaryCurvatureEnergy::BoundaryCurvatureEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetBoundaryCurvature>& targetBoundaryCurvatures)
    : ModularGaussNewtonEnergyBlock(targetBoundaryCurvatures.size(), ws.m_dofs.numDofs()), m_ws(ws), m_targetBoundaryCurvatures(targetBoundaryCurvatures)
{
    m_name = "BoundaryCurvatureEnergy";
}

void BoundaryCurvatureEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (const TargetBoundaryCurvature& target : m_targetBoundaryCurvatures)
    {
        const double k = computeBoundaryCurvature(m_ws, target.m_vertex);

        system.m_residual[ridx] = k - target.m_k;
        ridx++;
    }
}

void BoundaryCurvatureEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row = system.createJacobianRow();
    for (const TargetBoundaryCurvature& target : m_targetBoundaryCurvatures)
    {
        LocalGradient<Eigen::RowVector3d> gk;
        const double k = computeBoundaryCurvatureGradient(m_ws, target.m_vertex, gk);

        system.m_residual[ridx] = k - target.m_k;

        row.clearPattern();
        row.beginFill();
        for (pmp::Vertex vv : gk)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.insertSegment(dofIdx * 3, gk[vv]);
            }
        }
        row.endFill();
        system.insertJacobianRow(ridx, row);

        ridx++;
    }
}

void BoundaryCurvatureEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (const TargetBoundaryCurvature& target : m_targetBoundaryCurvatures)
    {
        weights[ridx] = energyWeight * m_ws.m_currentVertexDualAreas[target.m_vertex];
        ridx++;
    }
}
