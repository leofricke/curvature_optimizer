#include "PositionEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"

PositionEnergy::PositionEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetPosition>& targetPositions)
    : ModularGaussNewtonEnergyBlock(targetPositions.size() * 3, ws.m_dofs.numDofs()), m_ws(ws), m_targetPositions(targetPositions)
{
    m_name = "PositionEnergy";
}

void PositionEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (const TargetPosition& target : m_targetPositions)
    {
        const pmp::Vertex v0 = target.m_bary.aVertex(m_ws.m_mesh);
        const pmp::Vertex v1 = target.m_bary.bVertex(m_ws.m_mesh);
        const pmp::Vertex v2 = target.m_bary.cVertex(m_ws.m_mesh);
        const Eigen::Vector3d pos = target.m_bary.interpolate(m_ws.position(v0), m_ws.position(v1), m_ws.position(v2));
        system.m_residual.segment<3>(ridx) = pos - target.m_position;

        ridx += 3;
    }
}

void PositionEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    spmat::SpVec<double> row = system.createJacobianRow();
    for (const TargetPosition& target : m_targetPositions)
    {
        const pmp::Vertex v0 = target.m_bary.aVertex(m_ws.m_mesh);
        const pmp::Vertex v1 = target.m_bary.bVertex(m_ws.m_mesh);
        const pmp::Vertex v2 = target.m_bary.cVertex(m_ws.m_mesh);
        const Eigen::Vector3d pos = target.m_bary.interpolate(m_ws.position(v0), m_ws.position(v1), m_ws.position(v2));
        system.m_residual.segment<3>(ridx) = pos - target.m_position;

        const int dofIdx0 = m_ws.m_dofs.dofIdx(v0);
        const int dofIdx1 = m_ws.m_dofs.dofIdx(v1);
        const int dofIdx2 = m_ws.m_dofs.dofIdx(v2);
        for (size_t i = 0; i < 3; ++i)
        {
            row.clearPattern();
            row.beginFill();
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
            row.endFill();
            system.insertJacobianRow(ridx + i, row);
        }

        ridx += 3;
    }
}

void PositionEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (const TargetPosition& target : m_targetPositions)
    {
        weights.segment<3>(ridx).setConstant(energyWeight * m_ws.m_currentFaceAreas[target.m_bary.face(m_ws.m_mesh)] / std::pow(m_ws.m_currentMeshArea, 2.0));
        ridx += 3;
    }
}
