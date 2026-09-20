#include <gtest/gtest.h>

#include "curvature_optimizer/energies/CurvatureVariationEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class CurvatureVariationEnergyTest : public testing::Test
{
public:
    CurvatureVariationEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("torus_part.off"));

        ws.m_vertexCurvatureOperator = VertexCurvatureOperator::MeanGaussian;
        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(CurvatureVariationEnergyTest, TestCurvatureVariationEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
    ModularGaussNewtonEnergyFunction energyFunc({ std::make_shared<CurvatureVariationEnergy>(ws) });
    energyFunc.m_system.m_buildFullJacobian = true;
    energyFunc.evaluateResidualAndJacobian(dofs);

    const double eps = 1e-8;
    const double atol = 1e-4;
    for (size_t i = 0; i < energyFunc.m_numResiduals; ++i)
    {
        if (energyFunc.m_system.m_residual[i] <= 1e-5)
        {
            continue;
        }
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