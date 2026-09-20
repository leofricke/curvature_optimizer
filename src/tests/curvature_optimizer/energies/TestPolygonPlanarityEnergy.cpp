#include <gtest/gtest.h>

#include "curvature_optimizer/energies/PolygonPlanarityEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class PolygonPlanarityEnergyTest : public testing::Test
{
public:
    PolygonPlanarityEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("two_quads_bend.obj"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(PolygonPlanarityEnergyTest, TestPolygonPlanarityEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ModularGaussNewtonEnergyFunction energyFunc({ std::make_shared<PolygonPlanarityEnergy>(ws) });
    energyFunc.m_system.m_buildFullJacobian = true;
    energyFunc.evaluateResidualAndJacobian(dofs);

    const double eps = 1e-8;
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
            EXPECT_TRUE(checkDeriv(energyFunc.m_system.m_jacobian.coeff(i, j), d));
        }
    }
}
