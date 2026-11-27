#include "PolygonPlanarityEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../quantities/FaceAreaDerivatives.h"
#include "EnergyUtil.h"

PolygonPlanarityEnergy::PolygonPlanarityEnergy(const CurvatureOptimizerWorkspace& ws)
    : GaussNewtonEnergyModule("PolygonPlanarityEnergy", ws.m_dofs.numDofs()), m_ws(ws)
{
    for (pmp::Face f : m_ws.m_mesh.faces())
    {
        if (m_ws.m_mesh.valence(f) > 3)
        {
            m_polygons.push_back(f);
        }
    }
}

size_t PolygonPlanarityEnergy::numElements() const
{
    return m_polygons.size();
}

size_t PolygonPlanarityEnergy::numResiduals(size_t element) const
{
    return m_ws.m_mesh.valence(m_polygons[element]);
}

void PolygonPlanarityEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const pmp::Face f = m_polygons[element];
    const size_t val = m_ws.m_mesh.valence(f);

    LocalGradient<Eigen::Matrix3d> Jn;
    const Eigen::Vector3d n = ctx.withJacobian() ? computeFaceNormalJacobian(m_ws, f, Jn) : m_ws.m_currentFaceNormals[f];

    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
    {
        c += m_ws.position(v);
    }
    c /= static_cast<double>(val);

    size_t k = 0;
    for (pmp::Vertex v : m_ws.m_mesh.vertices(f))
    {
        const Eigen::Vector3d p = m_ws.position(v);
        const double r = (p - c).dot(n);
        ctx.residual(k) = std::sqrt(std::max(1e-10, std::abs(r)));

        if (ctx.withJacobian())
        {
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
            insertVertexGradient(ctx.jacobianRow(k), m_ws, g, c0);
        }

        k++;
    }
}

void PolygonPlanarityEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        residualWeights(element).setConstant(std::sqrt(m_ws.m_currentFaceAreas[m_polygons[element]]) / m_ws.m_currentMeshArea);
    }
}
