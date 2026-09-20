#include <gtest/gtest.h>

#include "curvature_optimizer/energies/CurvatureEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/TargetData.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class CurvatureEnergyTest : public testing::Test
{
public:
    CurvatureEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.m_vertexCurvatureOperator = VertexCurvatureOperator::MeanGaussian;
        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});

        for (pmp::Vertex v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                targetVertexCurvatures.emplace_back(v, 0.0, 0.0);
            }
        }
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
    std::vector<TargetVertexCurvature> targetVertexCurvatures;
};

TEST_F(CurvatureEnergyTest, TestCurvatureEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
    ModularGaussNewtonEnergyFunction energyFunc({ std::make_shared<CurvatureEnergy>(ws, targetVertexCurvatures) });
    energyFunc.m_system.m_buildFullJacobian = true;
    energyFunc.evaluateResidualAndJacobian(dofs);

    const double eps = 1e-7;
    const double atol = 1e-6;
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
