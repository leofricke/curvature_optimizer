#include <gtest/gtest.h>

#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/TargetData.h"
#include "curvature_optimizer/MeanGaussianVertexCurvature.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class MeanGaussianVertexCurvatureTest : public testing::Test
{
public:
    MeanGaussianVertexCurvatureTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("torus_part.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(MeanGaussianVertexCurvatureTest, ComputePrincipalCurvaturesGradient)
{
    const double eps = 1e-7;
    const double atol = 1e-6;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    LocalGradient<Eigen::RowVector3d, 8> gkmax;
    LocalGradient<Eigen::RowVector3d, 8> gkmin;
    for (pmp::Vertex v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
        {
            continue;
        }

        gkmax.clear();
        gkmin.clear();
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);

        MeanGaussianVertexCurvature curv;
        curv.computePrincipalCurvaturesGradient(ws, v, gkmax, gkmin);

        for (pmp::Vertex vv : gkmax)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const Eigen::Vector2d plus = curv.computePrincipalCurvatures(ws, v);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const Eigen::Vector2d minus = curv.computePrincipalCurvatures(ws, v);

                const Eigen::Vector2d d = (plus - minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(gkmax[vv][j], d[0], atol));
                EXPECT_TRUE(checkDeriv(gkmin[vv][j], d[1], atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
