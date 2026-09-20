#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/BoundaryCurvatureDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class BoundaryCurvatureDerivativesTest : public testing::Test
{
public:
    BoundaryCurvatureDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(BoundaryCurvatureDerivativesTest, ComputeBoundaryCurvatureGradient)
{
    const double eps = 1e-5;
    const double atol = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Vertex v0 : mesh.vertices())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d> g;
        computeBoundaryCurvatureGradient(ws, v0, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double K_plus = computeBoundaryCurvature(ws, v0);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double K_minus = computeBoundaryCurvature(ws, v0);

                const double dK = (K_plus - K_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dK, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
