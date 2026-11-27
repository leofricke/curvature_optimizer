#include <gtest/gtest.h>

#include "curvature_optimizer/energies/BoundaryCurvatureEnergy.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/TargetData.h"

#include "TestFilesUtil.h"
#include "DerivativeTestUtil.h"
#include <pmp/io/io.h>

class BoundaryCurvatureEnergyTest : public testing::Test
{
public:
    BoundaryCurvatureEnergyTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});

        for (pmp::Vertex v : mesh.vertices())
        {
            if (mesh.is_boundary(v))
            {
                targetBoundaryCurvatures.emplace_back(v, 0.0);
            }
        }
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
    std::vector<TargetBoundaryCurvature> targetBoundaryCurvatures;
};

TEST_F(BoundaryCurvatureEnergyTest, TestBoundaryCurvatureEnergyJacobian)
{
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
    GaussNewtonEnergyFunction energyFunc({ std::make_shared<BoundaryCurvatureEnergy>(ws, targetBoundaryCurvatures) });
    energyFunc.m_buildFullJacobian = true;
    energyFunc.evaluateEnergyAndGradient(dofs);

    const double eps = 1e-5;
    const double atol = 1e-7;
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
            EXPECT_TRUE(checkDeriv(energyFunc.m_jacobian.coeff(static_cast<int>(i), static_cast<int>(j)), d, atol));
        }
    }
}
