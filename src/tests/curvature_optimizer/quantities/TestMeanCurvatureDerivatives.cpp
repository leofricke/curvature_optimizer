#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/MeanCurvatureDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class MeanCurvatureDerivativesTest : public testing::Test
{
public:
    MeanCurvatureDerivativesTest()
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

TEST_F(MeanCurvatureDerivativesTest, ComputeVertexMeanCurvatureGradient)
{
    const double eps = 1e-8;
    const double atol = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Vertex v0 : mesh.vertices())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d, 8> g;
        computeVertexMeanCurvatureGradient(ws, v0, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double H_plus = computeVertexMeanCurvature(ws, v0);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double H_minus = computeVertexMeanCurvature(ws, v0);

                const double dH = (H_plus - H_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dH, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
