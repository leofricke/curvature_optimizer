#include "RegularizationEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"

RegularizationEnergy::RegularizationEnergy(const CurvatureOptimizerWorkspace& ws)
    : GaussNewtonEnergyModule("RegularizationEnergy", ws.m_dofs.numDofs()), m_ws(ws)
{
    resetPositions();
}

void RegularizationEnergy::resetPositions()
{
    m_initialPositions.clear();
    m_initialPositions.reserve(m_ws.m_mesh.n_vertices());
    for (pmp::Vertex v : m_ws.m_mesh.vertices())
    {
        m_initialPositions.emplace_back(v, m_ws.position(v));
    }
}

size_t RegularizationEnergy::numElements() const
{
    return m_initialPositions.size();
}

size_t RegularizationEnergy::numResiduals(size_t /*element*/) const
{
    return 3;
}

void RegularizationEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const InitialPosition& initialPosition = m_initialPositions[element];
    const Eigen::Vector3d p = m_ws.position(initialPosition.m_vertex);
    for (size_t i = 0; i < 3; ++i)
    {
        ctx.residual(i) = p[i] - initialPosition.m_position[i];
    }

    if (ctx.withJacobian())
    {
        const int dofIdx = m_ws.m_dofs.dofIdx(initialPosition.m_vertex);
        if (dofIdx >= 0)
        {
            for (int i = 0; i < 3; ++i)
            {
                ctx.jacobianRow(i).insert(dofIdx * 3 + i, 1.0);
            }
        }
    }
}

void RegularizationEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(m_ws.m_currentVertexDualAreas[m_initialPositions[element].m_vertex] / std::pow(m_ws.m_currentMeshArea, 2.0));
    }
}
