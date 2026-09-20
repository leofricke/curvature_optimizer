#include <gtest/gtest.h>

#include "curvature_optimizer/energies/ConformalEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/TargetData.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class ConformalEnergyTest : public testing::Test
{
public:
    ConformalEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});

        for (pmp::Face f : mesh.faces())
        {
            targetMetrics.push_back(createTargetMetric(ws, f));
            targetMetrics.back().m_triangles[0].m_le01 *= 2.0;
            targetMetrics.back().m_triangles[0].m_le12 *= 2.0;
            targetMetrics.back().m_triangles[0].m_le20 *= 2.0;
            targetMetrics.back().m_triangles[0].m_area *= 4.0;
        }
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
    std::vector<TargetMetric> targetMetrics;
};

TEST_F(ConformalEnergyTest, EquivalentToCotanEnergy)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ModularGaussNewtonEnergyFunction energyFunc({ std::make_shared<ConformalEnergy>(ws, targetMetrics) });
    energyFunc.evaluateEnergy(dofs);

    size_t idx = 0;
    for (const TargetMetric& target : targetMetrics)
    {
        const Eigen::Vector3d center = target.m_subdivision.getCenterPoint(ws, target.m_face);

        for (const VirtualTriangle& tri : target.m_triangles)
        {
            const double energy = energyFunc.m_system.m_residual.segment<6>(idx).squaredNorm();

            const pmp::Vertex v0 = ws.m_mesh.from_vertex(tri.m_halfedge);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(tri.m_halfedge);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);
            const Eigen::Vector3d e01 = p1 - p0;
            const Eigen::Vector3d e12 = center - p1;
            const Eigen::Vector3d e20 = p0 - center;
            const Eigen::Vector3d n = e01.cross(-e20);
            const double doubleArea = n.norm();

            const double cot0 = 0.5 * (std::pow(tri.m_le01, 2.0) - std::pow(tri.m_le12, 2.0) + std::pow(tri.m_le20, 2.0)) / (2.0 * tri.m_area);
            const double cot1 = 0.5 * (std::pow(tri.m_le01, 2.0) + std::pow(tri.m_le12, 2.0) - std::pow(tri.m_le20, 2.0)) / (2.0 * tri.m_area);
            const double cot2 = 0.5 * (-std::pow(tri.m_le01, 2.0) + std::pow(tri.m_le12, 2.0) + std::pow(tri.m_le20, 2.0)) / (2.0 * tri.m_area);
            const double cotanEnergy = (cot0 * e12.squaredNorm() + cot1 * e20.squaredNorm() + cot2 * e01.squaredNorm()) / doubleArea;

            EXPECT_NEAR(energy, cotanEnergy, 1e-10);

            idx += 6;
        }
    }
}

TEST_F(ConformalEnergyTest, TestConformalEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
    ModularGaussNewtonEnergyFunction energyFunc({ std::make_shared<ConformalEnergy>(ws, targetMetrics) });
    energyFunc.m_system.m_buildFullJacobian = true;
    energyFunc.evaluateResidualAndJacobian(dofs);

    const double eps = 1e-8;
    const double atol = 1e-7;
    for (size_t i = 0; i < energyFunc.m_numResiduals; ++i)
    {
        for (size_t j = 0; j < static_cast<size_t>(dofs.size()); ++j)
        {
            const double orig = dofs[j];

            dofs[j] = orig + eps;
            ws.updateGeometry();
            energyFunc.evaluateResidual(dofs);
            const Eigen::VectorXd r_plus = energyFunc.m_system.m_residual;

            dofs[j] = orig - eps;
            ws.updateGeometry();
            energyFunc.evaluateResidual(dofs);
            const Eigen::VectorXd r_minus = energyFunc.m_system.m_residual;

            dofs[j] = orig;
            ws.updateGeometry();

            const double d = (r_plus[i] - r_minus[i]) / (2.0 * eps);
            EXPECT_TRUE(checkDeriv(energyFunc.m_system.m_jacobian.coeff(i, j), d, atol));
        }
    }
}
