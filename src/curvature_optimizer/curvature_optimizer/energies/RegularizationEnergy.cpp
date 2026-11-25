#include "RegularizationEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"

RegularizationEnergy::RegularizationEnergy(const CurvatureOptimizerWorkspace& ws)
    : ModularGaussNewtonEnergyBlock(ws.m_mesh.n_vertices() * 3, ws.m_dofs.numDofs()), m_ws(ws)
{
    m_name = "RegularizationEnergy";

    m_initialPositions.reserve(m_ws.m_mesh.n_vertices());
    for (pmp::Vertex v : m_ws.m_mesh.vertices())
    {
        m_initialPositions.emplace_back(v, m_ws.position(v));
    }
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

void RegularizationEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (const InitialPosition& initialPosition : m_initialPositions)
    {
        const Eigen::Vector3d p = m_ws.position(initialPosition.m_vertex);
        system.m_residual.segment<3>(ridx) = p - initialPosition.m_position;

        ridx += 3;
    }
}

void RegularizationEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row = system.createJacobianRow();
    for (const InitialPosition& initialPosition : m_initialPositions)
    {
        const Eigen::Vector3d p = m_ws.position(initialPosition.m_vertex);
        system.m_residual.segment<3>(ridx) = p - initialPosition.m_position;

        const int dofIdx = m_ws.m_dofs.dofIdx(initialPosition.m_vertex);
        if (dofIdx >= 0)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                row.clearPattern();
                row.beginFill();
                row.insert(dofIdx * 3 + i, 1.0);
                row.endFill();
                system.insertJacobianRow(ridx + i, row);
            }
        }

        ridx += 3;
    }
}

void RegularizationEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (const InitialPosition& initialPosition : m_initialPositions)
    {
        weights.segment<3>(ridx).setConstant(energyWeight * m_ws.m_currentVertexDualAreas[initialPosition.m_vertex] / std::pow(m_ws.m_currentMeshArea, 2.0));
        ridx += 3;
    }
}
