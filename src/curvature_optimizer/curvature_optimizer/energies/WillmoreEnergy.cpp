#include "WillmoreEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/MeanCurvatureDerivatives.h"
#include "EnergyUtil.h"

WillmoreEnergy::WillmoreEnergy(CurvatureOptimizerWorkspace& ws)
    : GaussNewtonEnergyModule("WillmoreEnergy", ws.m_dofs.numDofs()), m_ws(ws)
{
}

size_t WillmoreEnergy::numElements() const
{
    return m_ws.m_mesh.vertices_size();
}

size_t WillmoreEnergy::numResiduals(size_t /*element*/) const
{
    return 1;
}

void WillmoreEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const pmp::Vertex v(static_cast<pmp::IndexType>(element));

    if (ctx.withJacobian())
    {
        LocalGradient<Eigen::RowVector3d, 8> gh;
        ctx.residual(0) = computeVertexMeanCurvatureGradient(m_ws, v, gh);
        insertVertexGradient(ctx.jacobianRow(0), m_ws, gh);
    }
    else
    {
        ctx.residual(0) = computeVertexMeanCurvature(m_ws, v);
    }
}

void WillmoreEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentVertexDualAreas[pmp::Vertex(static_cast<pmp::IndexType>(element))]);
    }
}
