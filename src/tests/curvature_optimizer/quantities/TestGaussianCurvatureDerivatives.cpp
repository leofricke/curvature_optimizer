#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/GaussianCurvatureDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class GaussianCurvatureDerivativesTest : public testing::Test
{
public:
    GaussianCurvatureDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(GaussianCurvatureDerivativesTest, ComputeVertexGaussianCurvatureGradient)
{
    const double eps = 1e-7;
    const double atol = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Vertex v : mesh.vertices())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d, 8> g;
        computeVertexGaussianCurvatureGradient(ws, v, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double K_plus = computeVertexGaussianCurvature(ws, v);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double K_minus = computeVertexGaussianCurvature(ws, v);

                const double dK = (K_plus - K_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dK, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
