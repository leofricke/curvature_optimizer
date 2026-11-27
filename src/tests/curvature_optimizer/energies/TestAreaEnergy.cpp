#include <gtest/gtest.h>

#include "curvature_optimizer/energies/AreaEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/TargetData.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class AreaEnergyTest : public testing::Test
{
public:
    AreaEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});

        for (pmp::Face f : mesh.faces())
        {
            targetMetrics.push_back(createTargetMetric(mesh, f));
        }
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
    std::vector<TargetMetric> targetMetrics;
};

TEST_F(AreaEnergyTest, TestAreaEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
    GaussNewtonEnergyFunction energyFunc({ std::make_shared<AreaEnergy>(ws, targetMetrics) });
    energyFunc.m_buildFullJacobian = true;
    energyFunc.evaluateEnergyAndGradient(dofs);

    const double eps = 1e-8;
    for (size_t j = 0; j < static_cast<size_t>(dofs.size()); ++j)
    {
        const double orig = dofs[j];

        dofs[j] = orig + eps;
        ws.updateGeometry();
        energyFunc.evaluateEnergy(dofs);
        const Eigen::VectorXd r_plus = energyFunc.m_residual;

        dofs[j] = orig - eps;
        ws.updateGeometry();
        energyFunc.evaluateEnergy(dofs);
        const Eigen::VectorXd r_minus = energyFunc.m_residual;

        dofs[j] = orig;
        ws.updateGeometry();

        for (size_t i = 0; i < energyFunc.m_numResiduals; ++i)
        {
            const double d = (r_plus[i] - r_minus[i]) / (2.0 * eps);
            EXPECT_TRUE(checkDeriv(energyFunc.m_jacobian.coeff(static_cast<int>(i), static_cast<int>(j)), d));
        }
    }
}
