#include "AreaEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "../quantities/DiffRules.h"

AreaEnergy::AreaEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetMetric>& targetMetrics)
    : ModularGaussNewtonEnergyBlock(targetMetrics.size() * 2, ws.m_dofs.numDofs()), m_ws(ws), m_targetMetrics(targetMetrics)
{
    m_name = "AreaEnergy";
}

void AreaEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (const TargetMetric& target : m_targetMetrics)
    {
        const double A = m_ws.m_currentFaceAreas[target.m_face];
        system.m_residual[ridx] = std::sqrt(target.m_area / A);
        system.m_residual[ridx + 1] = std::sqrt(A / target.m_area);

        ridx += 2;
    }
}

void AreaEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row = system.createJacobianRow();
    for (const TargetMetric& target : m_targetMetrics)
    {
        const double A = m_ws.m_currentFaceAreas[target.m_face];

        system.m_residual[ridx] = std::sqrt(target.m_area / A);
        system.m_residual[ridx + 1] = std::sqrt(A / target.m_area);

        LocalGradient<Eigen::RowVector3d> ga = m_ws.m_currentFaceAreaGradients[target.m_face];

        LocalGradient<Eigen::RowVector3d> gai;
        applyChainRule1(A, [](double x) { return -1.0 / (2.0 * x * std::sqrt(x)); }, ga, gai);
        applyChainRule1Inplace(A, [](double x) { return 1.0 / (2.0 * std::sqrt(x)); }, ga);

        row.clearPattern();
        row.beginFill();
        for (pmp::Vertex vv : gai)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.insertSegment(dofIdx * 3, gai[vv] * std::sqrt(target.m_area));
            }
        }
        row.endFill();
        system.insertJacobianRow(ridx, row);

        row.clearPattern();
        row.beginFill();
        for (pmp::Vertex vv : ga)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.insertSegment(dofIdx * 3, ga[vv] / std::sqrt(target.m_area));
            }
        }
        row.endFill();
        system.insertJacobianRow(ridx + 1, row);

        ridx += 2;
    }
}

void AreaEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (const TargetMetric& target : m_targetMetrics)
    {
        weights.segment<2>(ridx).setConstant(energyWeight * m_ws.m_currentFaceAreas[target.m_face] / m_ws.m_currentMeshArea);
        ridx += 2;
    }
}
