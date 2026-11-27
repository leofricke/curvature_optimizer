#include "BoundaryCurvatureEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "../quantities/BoundaryCurvatureDerivatives.h"
#include "EnergyUtil.h"

BoundaryCurvatureEnergy::BoundaryCurvatureEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetBoundaryCurvature>& targetBoundaryCurvatures)
    : GaussNewtonEnergyModule("BoundaryCurvatureEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_targetBoundaryCurvatures(targetBoundaryCurvatures)
{
}

size_t BoundaryCurvatureEnergy::numElements() const
{
    return m_targetBoundaryCurvatures.size();
}

size_t BoundaryCurvatureEnergy::numResiduals(size_t /*element*/) const
{
    return 1;
}

void BoundaryCurvatureEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const TargetBoundaryCurvature& target = m_targetBoundaryCurvatures[element];

    if (ctx.withJacobian())
    {
        LocalGradient<Eigen::RowVector3d> gk;
        const double k = computeBoundaryCurvatureGradient(m_ws, target.m_vertex, gk);
        ctx.residual(0) = k - target.m_k;
        insertVertexGradient(ctx.jacobianRow(0), m_ws, gk);
    }
    else
    {
        ctx.residual(0) = computeBoundaryCurvature(m_ws, target.m_vertex) - target.m_k;
    }
}

void BoundaryCurvatureEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentVertexDualAreas[m_targetBoundaryCurvatures[element].m_vertex]);
    }
}
