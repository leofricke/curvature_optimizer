#include "PolygonPlanarityEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/FaceAreaDerivatives.h"

PolygonPlanarityEnergy::PolygonPlanarityEnergy(const CurvatureOptimizerWorkspace& ws)
    : ModularGaussNewtonEnergyBlock(0, ws.m_dofs.numDofs()), m_ws(ws)
{
    m_name = "PolygonPlanarityEnergy";
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const size_t val = m_ws.m_mesh.valence(f);
        if (val > 3)
        {
            m_numResiduals += val;
        }
    }
}

void PolygonPlanarityEnergy::evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const size_t val = m_ws.m_mesh.valence(f);
        if (val > 3)
        {
            const Eigen::Vector3d& n = m_ws.m_currentFaceNormals[f];

            Eigen::Vector3d c = Eigen::Vector3d::Zero();
            for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
            {
                const Eigen::Vector3d p = m_ws.position(v);
                c += p;
            }
            c /= static_cast<double>(val);

            for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
            {
                const Eigen::Vector3d p = m_ws.position(v);
                const double r = (p - c).dot(n);
                system.m_residual[ridx] = std::sqrt(std::max(1e-10, std::abs(r)));

                ridx++;
            }
        }
    }
}

void PolygonPlanarityEnergy::evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system)
{
    size_t ridx = m_offset;
    SparseVector row;
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const size_t val = m_ws.m_mesh.valence(f);
        if (val > 3)
        {
            LocalGradient<Eigen::Matrix3d> Jn;
            const Eigen::Vector3d n = computeFaceNormalJacobian(m_ws, f, Jn);

            Eigen::Vector3d c = Eigen::Vector3d::Zero();
            for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
            {
                const Eigen::Vector3d p = m_ws.position(v);
                c += p;
            }
            c /= static_cast<double>(val);

            for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
            {
                const Eigen::Vector3d p = m_ws.position(v);
                const double r = (p - c).dot(n);
                system.m_residual[ridx] = std::sqrt(std::max(1e-10, std::abs(r)));

                LocalGradient<Eigen::RowVector3d> g;
                for (pmp::Vertex vv : Jn)
                {
                    g[vv] += Jn[vv].transpose() * (p - c);
                }
                g[v] += n;
                for (pmp::Vertex vv : m_ws.m_mesh.vertices(f))
                {
                    g[vv] -= (Eigen::Matrix3d::Identity() / val) * n;
                }

                const double sign = r >= 0.0 ? 1.0 : -1.0;
                const double c0 = (std::abs(r) > 1e-10) ? sign / (2.0 * std::sqrt(std::abs(r))) : 0.0;
                g *= c0;

                row.clear();
                for (pmp::Vertex vv : g)
                {
                    const int dofIdx = m_ws.m_dofs.dofIdx(vv);
                    if (dofIdx >= 0)
                    {
                        row.addSegment(dofIdx * 3, g[vv]);
                    }
                }
                system.insertJacobianRow(ridx, row);

                ridx++;
            }
        }
    }
}

void PolygonPlanarityEnergy::initWeights(double energyWeight, Eigen::VectorXd& weights)
{
    size_t ridx = m_offset;
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        const size_t val = m_ws.m_mesh.valence(f);
        if (val > 3)
        {
            weights.segment(ridx, val).setConstant(energyWeight * std::sqrt(m_ws.m_currentFaceAreas[f]) / m_ws.m_currentMeshArea);
            ridx += val;
        }
    }
}
