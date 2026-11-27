#include "ConformalEnergy.h"

#include "../CurvatureOptimizerWorkspace.h"
#include "../TargetData.h"
#include "../Tolerances.h"
#include "../quantities/DiffRules.h"
#include "../quantities/FaceAreaDerivatives.h"
#include "EnergyUtil.h"

ConformalEnergy::ConformalEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetMetric>& targetMetrics)
    : GaussNewtonEnergyModule("ConformalEnergy", ws.m_dofs.numDofs()), m_ws(ws), m_targetMetrics(targetMetrics)
{
}

size_t ConformalEnergy::numElements() const
{
    return m_targetMetrics.size();
}

size_t ConformalEnergy::numResiduals(size_t element) const
{
    return 6 * m_targetMetrics[element].m_triangles.size();
}

void ConformalEnergy::evaluateElement(size_t element, GaussNewtonElementContext& ctx) const
{
    const TargetMetric& target = m_targetMetrics[element];
    const Eigen::Vector3d center = target.m_subdivision.getCenterPoint(m_ws, target.m_face);

    size_t ridx = 0;
    for (const VirtualTriangle& tri : target.m_triangles)
    {
        const pmp::Vertex v0 = m_ws.m_mesh.from_vertex(tri.m_halfedge);
        const pmp::Vertex v1 = m_ws.m_mesh.to_vertex(tri.m_halfedge);
        const Eigen::Vector3d p0 = m_ws.position(v0);
        const Eigen::Vector3d p1 = m_ws.position(v1);
        const Eigen::Vector3d e01 = p1 - p0;
        const Eigen::Vector3d e12 = center - p1;
        const Eigen::Vector3d e20 = p0 - center;
        Eigen::Vector3d n = e01.cross(-e20);
        const double A = 0.5 * n.norm();
        n.normalize();
        const double a01 = std::sqrt(tri.m_area / (3.0 * std::pow(tri.m_le01, 2.0)));
        const double a12 = std::sqrt(tri.m_area / (3.0 * std::pow(tri.m_le12, 2.0)));
        const double a20 = std::sqrt(tri.m_area / (3.0 * std::pow(tri.m_le20, 2.0)));
        const double b = (tri.m_le01 + tri.m_le12 + tri.m_le20) / std::sqrt(12.0 * tri.m_area);

        const BarycentricCoords hb01 = tri.m_q01 - tri.m_incenter;
        const BarycentricCoords hb12 = tri.m_q12 - tri.m_incenter;
        const BarycentricCoords hb20 = tri.m_q20 - tri.m_incenter;
        const Eigen::Vector3d h01 = hb01.interpolate(p0, p1, center);
        const Eigen::Vector3d h12 = hb12.interpolate(p0, p1, center);
        const Eigen::Vector3d h20 = hb20.interpolate(p0, p1, center);

        // gradient of 1 / sqrt(A)
        LocalGradient<Eigen::RowVector3d> ga;
        if (ctx.withJacobian())
        {
            ga[v0] += 0.5 * n.cross(e12);
            ga[v1] += 0.5 * n.cross(e20);
            size_t hidx = 0;
            for (pmp::Halfedge h : m_ws.m_mesh.halfedges(target.m_face))
            {
                ga[m_ws.m_mesh.from_vertex(h)] += target.m_subdivision.m_centerWeights[hidx++] * 0.5 * n.cross(e01);
            }
            applyChainRule1Inplace(A, [](double x) { return -1.0 / (2.0 * x * std::sqrt(x)); }, ga);
        }

        // edge residuals
        for (size_t j = 0; j < 3; ++j)
        {
            const double aj = (j == 0) ? a01 : ((j == 1) ? a12 : a20);
            const Eigen::Vector3d& ej = (j == 0) ? e01 : ((j == 1) ? e12 : e20);
            const double l = ej.norm();

            ctx.residual(ridx + j) = aj * l / std::sqrt(A);

            if (ctx.withJacobian())
            {
                const Eigen::RowVector3d ej0 = ej.transpose().normalized();
                LocalGradient<Eigen::RowVector3d> gl;
                if (j == 0)
                {
                    gl[v1] += ej0;
                    gl[v0] += -ej0;
                }
                else if (j == 1)
                {
                    gl[v1] += -ej0;
                    size_t hidx = 0;
                    for (pmp::Halfedge h : m_ws.m_mesh.halfedges(target.m_face))
                    {
                        gl[m_ws.m_mesh.from_vertex(h)] += target.m_subdivision.m_centerWeights[hidx++] * ej0;
                    }
                }
                else
                {
                    gl[v0] += ej0;
                    size_t hidx = 0;
                    for (pmp::Halfedge h : m_ws.m_mesh.halfedges(target.m_face))
                    {
                        gl[m_ws.m_mesh.from_vertex(h)] += target.m_subdivision.m_centerWeights[hidx++] * -ej0;
                    }
                }

                LocalGradient<Eigen::RowVector3d> g;
                applyProductRule1(l, 1.0 / std::sqrt(A), gl, ga, g);
                insertVertexGradient(ctx.jacobianRow(ridx + j), m_ws, g, aj);
            }
        }

        // incircle residuals
        for (size_t j = 0; j < 3; ++j)
        {
            const Eigen::Vector3d& hj = (j == 0) ? h01 : ((j == 1) ? h12 : h20);
            const BarycentricCoords& hbj = (j == 0) ? hb01 : ((j == 1) ? hb12 : hb20);
            const double l = hj.norm();

            ctx.residual(ridx + 3 + j) = b * l / std::sqrt(A);

            if (ctx.withJacobian())
            {
                const Eigen::RowVector3d hj0 = hj.transpose().normalized();
                LocalGradient<Eigen::RowVector3d> gl;
                gl[v0] += hbj.m_alpha * hj0;
                gl[v1] += hbj.m_beta * hj0;
                size_t hidx = 0;
                for (pmp::Halfedge h : m_ws.m_mesh.halfedges(target.m_face))
                {
                    gl[m_ws.m_mesh.from_vertex(h)] += target.m_subdivision.m_centerWeights[hidx++] * hbj.m_gamma * hj0;
                }

                LocalGradient<Eigen::RowVector3d> g;
                applyProductRule1(l, 1.0 / std::sqrt(A), gl, ga, g);
                insertVertexGradient(ctx.jacobianRow(ridx + 3 + j), m_ws, g, b);
            }
        }

        if (n.dot(m_ws.m_currentFaceNormals[target.m_face]) < 0.0)
        {
            // flipped triangle: penalize the energy, but keep the gradient of the actual residuals
            if (ctx.withJacobian())
            {
                ctx.commitRows();
            }
            for (size_t k = 0; k < 6; ++k)
            {
                ctx.residual(ridx + k) = 1e7;
            }
        }

        ridx += 6;
    }
}

void ConformalEnergy::initialize()
{
    for (size_t element = 0; element < numElements(); ++element)
    {
        const TargetMetric& target = m_targetMetrics[element];
        auto weights = residualWeights(element);
        for (size_t t = 0; t < target.m_triangles.size(); ++t)
        {
            weights.segment<6>(static_cast<Eigen::Index>(6 * t)).setConstant(target.m_triangles[t].m_area / m_ws.m_currentMeshArea);
        }
    }
}
