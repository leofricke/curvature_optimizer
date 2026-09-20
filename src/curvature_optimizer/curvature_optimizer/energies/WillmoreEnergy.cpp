#include "WillmoreEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/MeanCurvatureDerivatives.h"

WillmoreEnergy::WillmoreEnergy(CurvatureOptimizerWorkspace& ws)
    : ModularGaussNewtonEnergyBlock(ws.m_mesh.n_vertices(), ws.m_dofs.numDofs()), m_ws(ws)
{
    m_name = "WillmoreEnergy";
}

void WillmoreEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (pmp::Vertex v : m_ws.m_mesh.vertices())
    {
        const double H = computeVertexMeanCurvature(m_ws, v);

        system.m_residual[ridx] = H;
        ridx++;
    }
}

void WillmoreEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    SparseVector row;
    for (pmp::Vertex v : m_ws.m_mesh.vertices())
    {
        LocalGradient<Eigen::RowVector3d, 8> gh;
        const double H = computeVertexMeanCurvatureGradient(m_ws, v, gh);

        system.m_residual[ridx] = H;

        row.clear();
        for (pmp::Vertex vv : gh)
        {
            const int dofIdx = m_ws.m_dofs.dofIdx(vv);
            if (dofIdx >= 0)
            {
                row.addSegment(dofIdx * 3, gh[vv]);
            }
        }
        system.insertJacobianRow(ridx, row);

        ridx++;
    }
}

void WillmoreEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (pmp::Vertex v : m_ws.m_mesh.vertices())
    {
        weights[ridx] = energyWeight * m_ws.m_currentVertexDualAreas[v];
        ridx++;
    }
}
