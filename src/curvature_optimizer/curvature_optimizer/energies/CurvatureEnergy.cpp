#include "CurvatureEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "EnergyUtil.h"

CurvatureEnergy::CurvatureEnergy(CurvatureOptimizerWorkspace& ws, const std::vector<TargetVertexCurvature>& targetVertexCurvatures)
    : GaussNewtonEnergyModule("CurvatureEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_targetVertexCurvatures(targetVertexCurvatures)
{
}

size_t CurvatureEnergy::numElements() const
{
    return m_targetVertexCurvatures.size();
}

size_t CurvatureEnergy::numResiduals(size_t /*element*/) const
{
    return 2;
}

void CurvatureEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const TargetVertexCurvature& target = m_targetVertexCurvatures[element];
    const PrincipalCurvatures& curv = m_ws.getCurrentVertexCurvature(target.m_vertex).m_curv;
    ctx.residual(0) = curv.m_kmax - target.m_kmax;
    ctx.residual(1) = curv.m_kmin - target.m_kmin;

    if (ctx.withJacobian())
    {
        insertVertexGradient(ctx.jacobianRow(0), m_ws, m_ws.m_currentMaxCurvatureGradients[target.m_vertex]);
        insertVertexGradient(ctx.jacobianRow(1), m_ws, m_ws.m_currentMinCurvatureGradients[target.m_vertex]);
    }
}

void CurvatureEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentVertexDualAreas[m_targetVertexCurvatures[element].m_vertex]);
    }
}
