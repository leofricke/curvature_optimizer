#include "AreaEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "../quantities/DiffRules.h"
#include "EnergyUtil.h"

AreaEnergy::AreaEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetMetric>& targetMetrics)
    : GaussNewtonEnergyModule("AreaEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_targetMetrics(targetMetrics)
{
}

size_t AreaEnergy::numElements() const
{
    return m_targetMetrics.size();
}

size_t AreaEnergy::numResiduals(size_t /*element*/) const
{
    return 2;
}

void AreaEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const TargetMetric& target = m_targetMetrics[element];
    const double A = m_ws.m_currentFaceAreas[target.m_face];

    ctx.residual(0) = std::sqrt(target.m_area / A);
    ctx.residual(1) = std::sqrt(A / target.m_area);

    if (ctx.withJacobian())
    {
        LocalGradient<Eigen::RowVector3d> ga = m_ws.m_currentFaceAreaGradients[target.m_face];

        LocalGradient<Eigen::RowVector3d> gai;
        applyChainRule1(A, [](double x) { return -1.0 / (2.0 * x * std::sqrt(x)); }, ga, gai);
        applyChainRule1Inplace(A, [](double x) { return 1.0 / (2.0 * std::sqrt(x)); }, ga);

        insertVertexGradient(ctx.jacobianRow(0), m_ws, gai, std::sqrt(target.m_area));
        insertVertexGradient(ctx.jacobianRow(1), m_ws, ga, 1.0 / std::sqrt(target.m_area));
    }
}

void AreaEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentFaceAreas[m_targetMetrics[element].m_face] / m_ws.m_currentMeshArea);
    }
}
