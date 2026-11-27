#include "PositionEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"

PositionEnergy::PositionEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetPosition>& targetPositions)
    : GaussNewtonEnergyModule("PositionEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_targetPositions(targetPositions)
{
}

size_t PositionEnergy::numElements() const
{
    return m_targetPositions.size();
}

size_t PositionEnergy::numResiduals(size_t /*element*/) const
{
    return 3;
}

void PositionEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const TargetPosition& target = m_targetPositions[element];
    const pmp::Vertex v0 = target.m_bary.aVertex(m_ws.m_mesh);
    const pmp::Vertex v1 = target.m_bary.bVertex(m_ws.m_mesh);
    const pmp::Vertex v2 = target.m_bary.cVertex(m_ws.m_mesh);
    const Eigen::Vector3d pos = target.m_bary.interpolate(m_ws.position(v0), m_ws.position(v1), m_ws.position(v2));
    for (size_t i = 0; i < 3; ++i)
    {
        ctx.residual(i) = pos[i] - target.m_position[i];
    }

    if (ctx.withJacobian())
    {
        const int dofIdx0 = m_ws.m_dofs.dofIdx(v0);
        const int dofIdx1 = m_ws.m_dofs.dofIdx(v1);
        const int dofIdx2 = m_ws.m_dofs.dofIdx(v2);
        for (int i = 0; i < 3; ++i)
        {
            JacobianRow& row = ctx.jacobianRow(i);
            if (dofIdx0 >= 0)
            {
                row.insert(dofIdx0 * 3 + i, target.m_bary.m_alpha);
            }
            if (dofIdx1 >= 0)
            {
                row.insert(dofIdx1 * 3 + i, target.m_bary.m_beta);
            }
            if (dofIdx2 >= 0)
            {
                row.insert(dofIdx2 * 3 + i, target.m_bary.m_gamma);
            }
        }
    }
}

void PositionEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentFaceAreas[m_targetPositions[element].m_bary.face(m_ws.m_mesh)] / std::pow(m_ws.m_currentMeshArea, 2.0));
    }
}
